#include <algorithm>
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <string>
#include <exception>
#include <cmath>
#include <cstdint>
#include <cstring>
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

std::size_t prec = 5;
#define BASE2_TYPE_SIZE 64
#define BASE2_TYPE_MASK 0x00000000FFFFFFFFUL

#define LOG_2_10 3.3219280948873623478703194294894

class fixed_point_exception : public std::exception{
private:
    const char * error_str;
public:
    fixed_point_exception(const char *str = nullptr) : error_str(str){}
    ~fixed_point_exception() throw(){}
    virtual const char *what() const throw(){
        if(error_str == nullptr){
            return "empty";
        }else{
            return error_str;
        }
    }
};

struct fixed_point{
    const std::size_t precision;
    int32_t sign;
    uint32_t *data;

    fixed_point(std::size_t prec, int n) : precision(prec), sign(n > 0 ? 1 : n < 0 ? -1 : 0), data(new uint32_t[prec]){
        for(size_t i = 0; i < prec; ++i){ data[i] = 0; }
        data[0] = static_cast<uint32_t>(n < 0 ? -n : n);
    }

    fixed_point(std::size_t prec, const std::string &str) : precision(prec), sign(1), data(new uint32_t[prec]){
        std::size_t n = static_cast<std::size_t>(std::ceil(LOG_2_10 * str.size() / (BASE2_TYPE_SIZE / 2)));
        if(n > prec){ throw(fixed_point_exception("input value is too large.")); }
        for(size_t i = 0; i < prec; ++i){ data[i] = 0; }
        fixed_point digit(prec, 1);
        char temp[] = { 0, '\0' };
        for(std::size_t i = 0; i < str.size(); ++i){
            primitive_mul_by_single(10);
            temp[0] = str[i];
            digit.data[0] = std::strtol(temp, nullptr, 10);
            primitive_add_n(digit);
        }
    }

    ~fixed_point(){
        delete[] data;
    }

    std::string to_string() const{
        if(sign == 0){ return "0"; }
        std::string result;
        fixed_point ten(precision, 10), temp(precision, 0);
        temp.sign = sign;
        for(std::size_t i = 0; i < precision; ++i){ temp.data[i] = data[i]; }
        for(std::size_t i = 0; !temp.primitive_check_zero(); ++i){
            result.push_back(static_cast<char>(temp.primitive_mod_by_single(10) + '0'));
            temp.primitive_div_by_single(10);
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

private:
    void primitive_add_n(const fixed_point &w){
        uint64_t carry2 = 0;
        for(std::size_t i = 0; i < precision; ++i){
            carry2 += static_cast<uint64_t>(data[i]) + w.data[i];
            data[i] = static_cast<uint64_t>(carry2 & BASE2_TYPE_MASK);
            carry2 >>= (BASE2_TYPE_SIZE / 2);
        }
    }

    void primitive_mul_by_single(uint32_t w){
        uint64_t carry2 = 0;
        for(std::size_t i = 0; i < precision; ++i){
            uint64_t x = static_cast<uint64_t>(data[i]) * w + carry2;
            data[i] = static_cast<uint32_t>(x & BASE2_TYPE_MASK);
            carry2 = x >> static_cast<uint64_t>(BASE2_TYPE_SIZE / 2);
        }
    }

    void primitive_div_by_single(uint32_t w){
        fixed_point result(precision, 1);
        uint64_t dividend = 0;
        for(int i = static_cast<int>(precision - 1); i >= 0; --i){
            dividend |= data[i];
            result.data[i] = static_cast<uint32_t>(dividend / w);
            dividend = (dividend % w) << (BASE2_TYPE_SIZE / 2);
        }
        delete data;
        data = result.data;
        result.data = nullptr;
    }

    uint32_t primitive_mod_by_single(uint32_t w) const{
        uint64_t p = 1, r = 0;
        for(std::size_t i = 0; i < precision; ++i){
            for(std::size_t bit = 0; bit < BASE2_TYPE_SIZE / 2; ++bit){
                if((data[i] & (1 << bit)) != 0){
                    r += p;
                    if(r >= w){ r -= w; }
                }
                p <<= 1;
                if(p >= w){ p -= w; }
            }
        }
        return static_cast<uint32_t>(r);
    }

    bool primitive_check_zero() const{
        for(std::size_t i = 0; i < precision; ++i){
            if(data[i] != 0){ return false; }
        }
        return true;
    }
};

int main(){
    cl_int err = CL_SUCCESS;
    try{
        std::ifstream ifile("cl/cl.cl", std::ios::binary | std::ios::ate);
        std::vector<char> cl_source;
        if(ifile.fail()){
            throw(std::exception("can not open the file \"cl/cl.cl\"."));
        }
        cl_source.resize((std::size_t)ifile.tellg());
        ifile.seekg(0, ifile.beg);
        ifile.read(&cl_source[0], cl_source.size());

        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if(platforms.size() == 0){
            std::cerr << "Platform size 0" << std::endl;
            return -1;
        }

        cl::Platform platform = platforms[0];
        cl_context_properties properties[] = {
            CL_CONTEXT_PLATFORM,
            reinterpret_cast<cl_context_properties>(platform()),
            0
        };
        cl::Context context(CL_DEVICE_TYPE_GPU, properties);

        std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
        cl::Device device = devices[0];
        cl::CommandQueue queue(context, device);

        cl::Program::Sources sources;
        sources.push_back(std::make_pair(&cl_source[0], cl_source.size()));
        cl::Program program(context, sources);
        program.build(devices);

        fixed_point f(prec, "123456789"), g(prec, "987654321"), h(prec, 1);
        std::cout << f.to_string() << std::endl;
        std::cout << g.to_string() << std::endl;

        cl::Buffer
            f_sign_buffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &f.sign),
            f_buffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(uint32_t) * prec, f.data),
            g_sign_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &g.sign),
            g_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(uint32_t) * prec, g.data);
        queue.enqueueWriteBuffer(f_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &f.sign);
        queue.enqueueWriteBuffer(f_buffer, CL_TRUE, 0, sizeof(uint32_t) * prec, f.data);
        queue.enqueueWriteBuffer(g_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &g.sign);
        queue.enqueueWriteBuffer(g_buffer, CL_TRUE, 0, sizeof(uint32_t) * prec, g.data);

        cl::Kernel kernel(program, "cl_main");
        kernel.setArg(0, f_sign_buffer);
        kernel.setArg(1, f_buffer);
        kernel.setArg(2, g_sign_buffer);
        kernel.setArg(3, g_buffer);

        cl::Event event;
        queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(1), cl::NDRange(1, 1), nullptr, &event);
        event.wait();
        queue.enqueueReadBuffer(f_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &f.sign, nullptr, &event);
        queue.enqueueReadBuffer(f_buffer, CL_TRUE, 0, sizeof(uint32_t) * prec, f.data, nullptr, &event);

        std::cout << f.to_string() << std::endl;
    }catch(cl::Error err){
        std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    }catch(std::exception err){
        std::cerr << "ERROR: " << err.what() << std::endl;
    }
    return 0;
}
