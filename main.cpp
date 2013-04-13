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

std::size_t g_integral_part = 1, g_fraction_part = 4, g_prec = g_integral_part + g_fraction_part;
#define BASE2_TYPE_SIZE 64
#define BASE2_TYPE_MASK 0x00000000FFFFFFFFUL

#define LOG_2_10 3.3219280948873623478703194294894

class fixed_point_exception : public std::exception{
private:
    const char *error_str;

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
    const std::size_t integral_part, fraction_part, precision;
    int32_t sign;
    uint32_t *data;

    fixed_point(std::size_t integral_part_, std::size_t fraction_part_, int n) :
        integral_part(integral_part_),
        fraction_part(fraction_part_),
        precision(integral_part_ + fraction_part_),
        sign(n > 0 ? 1 : n < 0 ? -1 : 0),
        data(new uint32_t[precision])
    {
        for(size_t i = 0; i < precision; ++i){ data[i] = 0; }
        data[fraction_part] = static_cast<uint32_t>(n < 0 ? -n : n);
    }

    fixed_point(
        std::size_t integral_part_,
        std::size_t fraction_part_,
        const std::string &str_i,
        const std::string &str_f
    ) :
        integral_part(integral_part_),
        fraction_part(fraction_part_),
        precision(integral_part_ + fraction_part_),
        sign(1),
        data(new uint32_t[precision])
    {
        for(std::size_t i = 0; i < precision; ++i){ data[i] = 0; }
        {
            // interpret integral part.
            std::size_t n = static_cast<std::size_t>(std::ceil(LOG_2_10 * str_i.size() / (BASE2_TYPE_SIZE / 2)));
            if(n > integral_part){ throw(fixed_point_exception("input value is too large.")); }
            fixed_point digit(integral_part, fraction_part, 1), z(integral_part, fraction_part, 1);
            for(size_t i = 0; i < precision; ++i){ z.data[i] = 0; }
            char temp[] = { 0, '\0' };
            for(std::size_t i = 0; i < str_i.size(); ++i){
                z.primitive_mul_by_single(10);
                temp[0] = str_i[i];
                digit.data[fraction_part] = std::strtol(temp, nullptr, 10);
                z.primitive_add_n(digit);
            }
            primitive_add_n(z);
        }
        {
            // interpret fraction part.
            fixed_point digit(integral_part, fraction_part, 1), z(integral_part, fraction_part, 1);
            for(size_t i = 0; i < precision; ++i){ z.data[i] = 0; }
            digit.data[fraction_part] = 0;
            char temp[] = { 0, '\0' };
            for(int i = static_cast<int>(str_f.size() - 1); i >= 0; --i){
                temp[0] = str_f[i];
                digit.data[fraction_part] = std::strtol(temp, nullptr, 10);
                z.primitive_add_n(digit);
                z.primitive_div_by_single(10);
            }
            primitive_add_n(z);
        }
    }

    ~fixed_point(){
        delete[] data;
    }

    std::string to_string() const{
        if(sign == 0){ return "0"; }
        std::string result;
        fixed_point temp(integral_part, fraction_part, 0);
        temp.sign = sign;
        if(sign < 0){ result.append("-"); }
        for(std::size_t i = 0; i < precision; ++i){ temp.data[i] = data[i]; }
        for(std::size_t i = 0; !temp.primitive_check_zero(); ++i){
            result.push_back(static_cast<char>(temp.primitive_mod_by_single(10) + '0'));
            temp.primitive_div_by_single(10);
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    std::string to_fp_string() const{
        if(sign == 0){ return "0"; }
        std::string result;
        {
            fixed_point integer(integral_part, 0, 0);
            integer.sign = 1;
            for(std::size_t i = 0; i < integral_part; ++i){
                integer.data[i] = data[fraction_part + i];
            }
            if(integer.primitive_check_zero()){
                result += "0";
            }else{
                result += integer.to_string();
            }
        }
        fixed_point temp(integral_part, fraction_part, 0);
        temp.sign = sign;
        for(std::size_t i = 0; i < fraction_part; ++i){ temp.data[i] = data[i]; }
        std::string str_f;
        if(temp.primitive_check_zero()){
            str_f = "0";
        }else{
            while(!temp.primitive_check_zero()){
                temp.primitive_mul_by_single(10);
                str_f.push_back(static_cast<char>(temp.data[fraction_part] + '0'));
                temp.data[fraction_part] = 0;
            }
        }
        result.push_back('.');
        result += str_f;
        return result;
    }

    static void fp_mul(fixed_point *z, const fixed_point *v, const fixed_point *w){
        if(v->sign == 0 || w->sign == 0){
            z->sign = 0;
            return;
        }
        z->sign = v->sign * w->sign;
        std::unique_ptr<uint32_t[]> buff_instance(new uint32_t[v->precision + w->precision]);
        uint32_t *buff = buff_instance.get();
        for(std::size_t i = 0; i < z->precision * 2; ++i){ buff[i] = 0; }
        for(std::size_t i = 0; i < v->precision; ++i){
            uint64_t carry2 = 0;
            for(std::size_t j = 0; j < w->precision; ++j){
                std::size_t k = i + j;
                uint64_t n = (uint64_t)(v->data[i]) * (uint64_t)(w->data[j]) + carry2;
                buff[k] += (uint32_t)(n & BASE2_TYPE_MASK);
                carry2 = (uint32_t)(n >> (BASE2_TYPE_SIZE / 2));
            }
            buff[i + w->precision] += (uint32_t)carry2;
        }
        for(std::size_t i = 0; i < z->precision; ++i){
            z->data[i] = buff[z->precision + i - 1];
        }
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
        fixed_point result(integral_part, fraction_part, 1);
        uint64_t dividend = 0;
        for(int i = static_cast<int>(precision - 1); i >= 0; --i){
            dividend |= data[i];
            result.data[i] = static_cast<uint32_t>(dividend / w);
            dividend = (dividend % w) << (BASE2_TYPE_SIZE / 2);
        }
        delete[] data;
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

        fixed_point f(g_integral_part, g_fraction_part, "1024201", "9876543210123456789");
        std::cout << f.to_fp_string() << std::endl;

        //fixed_point
        //    f(g_integral_part, g_fraction_part, 1),
        //    g(g_integral_part, g_fraction_part, 2),
        //    h(g_integral_part, g_fraction_part, 1);
        //f.data[4] = 1;
        //f.data[3] = 0x80000000;
        //f.data[0] = f.data[1] = f.data[2] = 0;
        //std::cout << f.to_fp_string() << std::endl;
        //std::cout << g.to_fp_string() << std::endl;

        //cl::Buffer
        //    f_sign_buffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &f.sign),
        //    f_buffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(uint32_t) * f.precision, f.data),
        //    g_sign_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &g.sign),
        //    g_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(uint32_t) * g.precision, g.data);
        //queue.enqueueWriteBuffer(f_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &f.sign);
        //queue.enqueueWriteBuffer(f_buffer, CL_TRUE, 0, sizeof(uint32_t) * f.precision, f.data);
        //queue.enqueueWriteBuffer(g_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &g.sign);
        //queue.enqueueWriteBuffer(g_buffer, CL_TRUE, 0, sizeof(uint32_t) * g.precision, g.data);

        //cl::Kernel kernel(program, "cl_main");
        //kernel.setArg(0, f_sign_buffer);
        //kernel.setArg(1, f_buffer);
        //kernel.setArg(2, g_sign_buffer);
        //kernel.setArg(3, g_buffer);

        //cl::Event event;
        //queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(1), cl::NDRange(1, 1), nullptr, &event);
        //event.wait();
        //queue.enqueueReadBuffer(f_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &f.sign, nullptr, &event);
        //queue.enqueueReadBuffer(f_buffer, CL_TRUE, 0, sizeof(uint32_t) * f.precision, f.data, nullptr, &event);

        //std::cout << f.to_fp_string() << std::endl;
    }catch(cl::Error err){
        std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    }catch(std::exception err){
        std::cerr << "ERROR: " << err.what() << std::endl;
    }
    return 0;
}
