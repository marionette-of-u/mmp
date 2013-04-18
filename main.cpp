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

class fixed_point{
private:
    std::size_t integral_part, fraction_part, precision;

public:
    int32_t sign;
    uint32_t *data;

public:
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

    std::string to_string(uint32_t radix = 10) const{
        if(sign == 0){ return "0"; }
        std::string result;
        fixed_point temp(integral_part, fraction_part, 0);
        temp.sign = sign;
        for(std::size_t i = 0; i < precision; ++i){ temp.data[i] = data[i]; }
        for(std::size_t i = 0; !temp.primitive_check_zero(); ++i){
            uint32_t n = temp.primitive_mod_by_single(radix);
            char c = 0;
            if(n >= 0 && n <= 9){
                c = static_cast<char>(n + '0');
            }else{
                c = static_cast<char>(n - 10 + 'A');
            }
            result.push_back(c);
            temp.primitive_div_by_single(radix);
        }
        if(sign < 0){ result.push_back('-'); }
        std::reverse(result.begin(), result.end());
        return result;
    }

    std::string to_fp_string(uint32_t radix = 10) const{
        if(sign == 0){ return "0"; }
        std::string result;
        if(sign < 0){
            result += "-";
        }
        {
            fixed_point integer(integral_part, 0, 0);
            integer.sign = 1;
            for(std::size_t i = 0; i < integral_part; ++i){
                integer.data[i] = data[fraction_part + i];
            }
            if(integer.primitive_check_zero()){
                result += "0";
            }else{
                result += integer.to_string(radix);
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
                temp.primitive_mul_by_single(radix);
                uint32_t n = temp.data[fraction_part];
                if(n >= 0 && n <= 9){
                    str_f.push_back(static_cast<char>(n + '0'));
                }else{
                    str_f.push_back(static_cast<char>(n - 10 + 'A'));
                }
                temp.data[fraction_part] = 0;
            }
        }
        result.push_back('.');
        result += str_f;
        return result;
    }

    void set(const fixed_point &w){
        primitive_set_fixed_point(w);
    }

    void resize(std::size_t i_part, std::size_t f_part){
        uint32_t *new_data = new uint32_t[i_part + f_part];
        if(f_part > fraction_part){
            for(std::size_t i = 0; i < fraction_part; ++i){
                new_data[f_part - fraction_part + i] = data[i];
            }
            for(std::size_t i = 0; i < f_part - fraction_part; ++i){
                new_data[i] = 0;
            }
        }else if(f_part < fraction_part){
            for(std::size_t i = 0; i < f_part; ++i){
                new_data[i] = data[fraction_part - f_part + i];
            }
        }else{
            for(std::size_t i = 0; i < fraction_part; ++i){
                new_data[i] = data[i];
            }
        }
        if(i_part > integral_part){
            for(std::size_t i = 0; i < integral_part; ++i){
                new_data[f_part + i] = data[fraction_part + i];
            }
            for(std::size_t i = 0; i < i_part - integral_part; ++i){
                new_data[f_part + integral_part + i] = 0;
            }
        }else if(i_part < integral_part){
            for(std::size_t i = 0; i < i_part; ++i){
                new_data[f_part + i] = data[fraction_part + i];
            }
        }else{
            for(std::size_t i = 0; i < integral_part; ++i){
                new_data[f_part + i] = data[fraction_part + i];
            }
        }
        fraction_part = f_part;
        integral_part = i_part;
        precision = fraction_part + integral_part;
        delete[] data;
        data = new_data;
    }

    void mul(const fixed_point &v, const fixed_point &w){
        if(v.sign == 0 || w.sign == 0){
            sign = 0;
            return;
        }
        sign = v.sign * w.sign;
        std::unique_ptr<uint32_t[]> buff_instance(new uint32_t[precision * 2]);
        uint32_t *buff = buff_instance.get();
        for(std::size_t i = 0; i < precision * 2; ++i){ buff[i] = 0; }
        for(std::size_t i = 0; i < v.precision; ++i){
            uint64_t carry2 = 0;
            for(std::size_t j = 0; j < w.precision; ++j){
                std::size_t k = i + j;
                uint64_t n = static_cast<uint64_t>(v.data[i]) * static_cast<uint64_t>(w.data[j]) + carry2;
                buff[k] += static_cast<uint32_t>(n & BASE2_TYPE_MASK);
                carry2 = static_cast<uint32_t>(n >> (BASE2_TYPE_SIZE / 2));
            }
            buff[i + w.precision] += static_cast<uint32_t>(carry2);
        }
        for(std::size_t i = 0; i < precision; ++i){
            data[i] = buff[precision + i - 1];
        }
    }

    void div(const fixed_point &u, const fixed_point &v){
        if(u.sign == 0){
            primitive_set_zero();
            return;
        }
        sign = u.sign * v.sign;
        const int32_t m = u.primitive_degree(), n = v.primitive_degree();
        const uint64_t b = static_cast<uint64_t>(BASE2_TYPE_MASK) + 1;
        uint32_t *un, *vn;
        uint64_t qhat, rhat, p;
        int64_t i, j, t, k;
        for(i = 0; i < static_cast<int64_t>(precision); ++i){ data[i] = 0; }
        const uint32_t s = primitive_nlz(v.data[n - 1]);
        std::unique_ptr<uint32_t[]> vn_scoped_guard(new uint32_t[precision * 2]);
        vn = vn_scoped_guard.get();
        for(i = n - 1; i > 0; --i){
            vn[i] = static_cast<uint32_t>(primitive_s_lshift(v.data[i], s) | primitive_s_rshift(v.data[i - 1], (BASE2_TYPE_SIZE / 2 - s)));
        }
        vn[0] = static_cast<uint32_t>(primitive_s_lshift(v.data[0], s));
        std::unique_ptr<uint32_t[]> un_scoped_guard(new uint32_t[(precision + 1) * 2 + precision]);
        un = un_scoped_guard.get();
        un[m + fraction_part] = static_cast<uint32_t>(primitive_s_rshift(u.data[m - 1], (BASE2_TYPE_SIZE / 2 - s)));
        for(i = 0; i < static_cast<int64_t>(fraction_part); ++i){ un[i] = 0; }
        for(i = m - 1; i > 0; --i){
            un[i + fraction_part] = static_cast<uint32_t>(primitive_s_lshift(u.data[i], s) | primitive_s_rshift(u.data[i - 1], (BASE2_TYPE_SIZE / 2 - s)));
        }
        un[fraction_part] = static_cast<uint32_t>(primitive_s_lshift(u.data[0], s));
        for(j = m + fraction_part - n; j >= 0; --j){
            qhat = (un[j + n] * b + un[j + n - 1]) / vn[n - 1];
            rhat = (un[j + n] * b + un[j + n - 1]) - qhat * vn[n - 1];
            do{
                if(qhat >= b || qhat * vn[n - 2] > b * rhat + un[j - n - 2]){
                    --qhat;
                    rhat += vn[n - 1];
                    if(rhat < b){ continue; }
                }
            }while(false);
            k = 0;
            for(i = 0; i < n; ++i){
                p = qhat * vn[i];
                t = un[i + j] - k - (p & BASE2_TYPE_MASK);
                un[i + j] = static_cast<uint32_t>(t);
                k = (p >> BASE2_TYPE_SIZE / 2) - (t >> BASE2_TYPE_SIZE / 2);
            }
            t = un[j + n] - k;
            un[j + n] = static_cast<uint32_t>(t);
            data[j] = static_cast<uint32_t>(qhat);
            if(t < 0){
                --data[j];
                k = 0;
                for(i = 0; i < n; ++i){
                    t = un[i + j] + vn[i] + k;
                    un[i + j] = static_cast<uint32_t>(t);
                }
                k = t >> BASE2_TYPE_SIZE / 2;
            }
        }
    }

    void add(const fixed_point &w){
        if(w.sign == 0){ return; }
        if(sign == 0){
            primitive_set_fixed_point(w);
            return;
        }
        if(sign == w.sign){
            primitive_add_n(w);
            return;
        }
        int u = primitive_compare(*this, w);
        if(u != 0){
            if(u > 0){
                primitive_sub_n(w);
            }else{
                fixed_point temp(integral_part, fraction_part, 1);
                temp.primitive_set_fixed_point(*this);
                primitive_set_fixed_point(w);
                primitive_sub_n(temp);
                sign = -sign;
            }
            primitive_zero_normalize();
        }else{
            primitive_set_zero();
        }
    }

    void sub(const fixed_point &w){
        if(w.sign == 0){ return; }
        if(sign == 0){
            primitive_set_fixed_point(w);
            sign = -sign;
        }
        if(sign != w.sign){
            primitive_add_n(w);
            return;
        }
        int u = primitive_compare(*this, w);
        if(u != 0){
            if(u > 0){
                primitive_sub_n(w);
            }else{
                fixed_point temp(integral_part, fraction_part, 1);
                temp.primitive_set_fixed_point(*this);
                primitive_set_fixed_point(w);
                primitive_sub_n(temp);
                sign = -sign;
            }
            primitive_zero_normalize();
        }else{
            primitive_set_zero();
        }
    }

private:
    void primitive_mul_inv(){
        if((data[0] & 1) == 0){ data[0] += 1; }
        fixed_point
            xn(integral_part, fraction_part, 0),
            t(integral_part, fraction_part, 0),
            u(integral_part, fraction_part, 0);
        xn.primitive_set_fixed_point(*this);
        while(true){
            t.primitive_mul_with_sign(*this, xn);
            std::cout << t.to_string() << std::endl;
            if(t.sign == 1){
                if(t.data[0] == 1){
                    bool f = true;
                    for(std::size_t i = 1; f && i < precision; ++i){
                        f = f && t.data[i] == 0;
                    }
                    if(f){
                        primitive_set_fixed_point(xn);
                        return;
                    }
                }
            }else if(t.sign == -1){
                bool f = true;
                for(std::size_t i = 0; f && i < precision; ++i){
                    f = f && t.data[i] == static_cast<uint32_t>(BASE2_TYPE_MASK);
                }
                if(f){
                    primitive_set_fixed_point(xn);
                    return;
                }
            }
            t.sign = -t.sign;
            t.primitive_add_n_by_single(2);
            u.primitive_mul_with_sign(xn, t);
            xn.primitive_set_fixed_point(u);
        }
    }

    void primitive_mul_with_sign(const fixed_point &v, const fixed_point &w){
        if(v.sign == 0 || w.sign == 0){
            sign = 0;
            return;
        }
        for(std::size_t i = 0; i < precision; ++i){ data[i] = 0; }
        sign = v.sign * w.sign;
        for(std::size_t i = 0; i < precision; ++i){
            uint64_t carry2 = 0;
            for(std::size_t j = 0; j < precision; ++j){
                std::size_t k = i + j;
                if(k >= precision){ continue; }
                uint64_t n = static_cast<uint64_t>(v.data[i]) * static_cast<uint64_t>(w.data[j]) + carry2;
                data[k] += static_cast<uint32_t>(n & BASE2_TYPE_MASK);
                carry2 = static_cast<uint32_t>(n >> (BASE2_TYPE_SIZE / 2));
            }
        }
    }

    int32_t primitive_degree() const{
        if(sign == 0){ return 0; }
        int32_t n = static_cast<int32_t>(precision - 1);
        while(data[n] == 0){ --n; }
        return n + 1;
    }

    void primitive_add_n(const fixed_point &w){
        uint64_t carry2 = 0;
        for(std::size_t i = 0; i < precision; ++i){
            carry2 += static_cast<uint64_t>(data[i]) + w.data[i];
            data[i] = static_cast<uint32_t>(carry2 & BASE2_TYPE_MASK);
            carry2 >>= (BASE2_TYPE_SIZE / 2);
        }
    }

    void primitive_add_n_by_single(uint32_t w){
        uint64_t carry2 = w;
        for(std::size_t i = 0; i < precision; ++i){
            carry2 += static_cast<uint64_t>(data[i]);
            data[i] = static_cast<uint32_t>(carry2 & BASE2_TYPE_MASK);
            carry2 >>= (BASE2_TYPE_SIZE / 2);
        }
    }

    void primitive_sub_n(const fixed_point &w){
        uint64_t borrow2 = 0;
        for(std::size_t i = 0; i < precision; ++i){
            uint64_t x = static_cast<uint64_t>(data[i]) - static_cast<uint64_t>(w.data[i]) - borrow2;
            data[i] = static_cast<uint32_t>(x & BASE2_TYPE_MASK);
            borrow2 = x > BASE2_TYPE_MASK ? 1 : 0;
        }
    }

    void primitive_sub_n_by_single(uint32_t w){
        uint64_t borrow2 = w;
        for(std::size_t i = 0; i < precision; ++i){
            uint64_t x = static_cast<uint64_t>(data[i]) - borrow2;
            data[i] = static_cast<uint32_t>(x & BASE2_TYPE_MASK);
            borrow2 = x > BASE2_TYPE_MASK ? 1 : 0;
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

    void primitive_set_fixed_point(const fixed_point &w){
        if(this == &w){ return; }
        sign = w.sign;
        for(std::size_t i = 0; i < precision; ++i){ data[i] = w.data[i]; }
    }

    int primitive_compare(const fixed_point &z, const fixed_point &w){
        for(int i = static_cast<int>(z.precision - 1); i >= 0; --i){
            if(z.data[i] > w.data[i]){ return 1; }
            if(z.data[i] < w.data[i]){ return -1; }
        }
        return 0;
    }

    void primitive_set_zero(){
        sign = 0;
    }

    void primitive_zero_normalize(){
        if(sign == 0){ return; }
        if(primitive_check_zero()){ sign = 0; }
    }

    static int32_t primitive_nlz(uint32_t x){
        uint32_t y;
        int32_t n = BASE2_TYPE_SIZE / 2, c = n / 2;
        do{
            y = x >> c;
            if(y != 0){
                n -= c;
                x = y;
            }
            c >>= 1;
        }while(c != 0);
        return n - static_cast<int32_t>(x);
    }

    static uint64_t primitive_s_lshift(uint64_t x, int32_t n){
        if(n >= BASE2_TYPE_SIZE / 2){ return 0; }
        return x << n;
    }

    static uint64_t primitive_s_rshift(uint64_t x, int32_t n){
        if(n >= BASE2_TYPE_SIZE / 2){ return 0; }
        return x >> n;
    }
};

int main(){
    fixed_point
        f(g_integral_part, g_fraction_part * 2, "1", "0"),
        g(g_integral_part, g_fraction_part * 2, "3", "0"),
        h(g_integral_part, g_fraction_part * 2, 0);
    std::cout << f.to_fp_string() << std::endl;
    std::cout << g.to_fp_string() << std::endl;

    // 0.333... = 0.5555555555555555555555555555555555555555555555555555555555555555
    h.div(f, g);
    std::cout << h.to_fp_string(0x10) << std::endl;

    // 0.999... = 0.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    f.set(h);
    h.mul(f, g);
    std::cout << h.to_fp_string(0x10) << std::endl;

    return 0;

    //cl_int err = CL_SUCCESS;
    //try{
    //    std::ifstream ifile("cl/cl.cl", std::ios::binary | std::ios::ate);
    //    std::vector<char> cl_source;
    //    if(ifile.fail()){
    //        throw(std::exception("can not open the file \"cl/cl.cl\"."));
    //    }
    //    cl_source.resize((std::size_t)ifile.tellg());
    //    ifile.seekg(0, ifile.beg);
    //    ifile.read(&cl_source[0], cl_source.size());

    //    std::vector<cl::Platform> platforms;
    //    cl::Platform::get(&platforms);
    //    if(platforms.size() == 0){
    //        std::cerr << "Platform size 0" << std::endl;
    //        return -1;
    //    }

    //    cl::Platform platform = platforms[0];
    //    cl_context_properties properties[] = {
    //        CL_CONTEXT_PLATFORM,
    //        reinterpret_cast<cl_context_properties>(platform()),
    //        0
    //    };
    //    cl::Context context(CL_DEVICE_TYPE_GPU, properties);

    //    std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
    //    cl::Device device = devices[0];
    //    cl::CommandQueue queue(context, device);

    //    cl::Program::Sources sources;
    //    sources.push_back(std::make_pair(&cl_source[0], cl_source.size()));
    //    cl::Program program(context, sources);
    //    program.build(devices, "-cl-strict-aliasing");

    //    fixed_point
    //        f(g_integral_part, g_fraction_part, "1", "0"),
    //        g(g_integral_part, g_fraction_part, "3", "0");
    //    std::cout << f.to_fp_string() << std::endl;
    //    std::cout << g.to_fp_string() << std::endl;

    //    cl::Buffer
    //        f_sign_buffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &f.sign),
    //        f_buffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(uint32_t) * f.precision, f.data),
    //        g_sign_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &g.sign),
    //        g_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(uint32_t) * g.precision, g.data);
    //    queue.enqueueWriteBuffer(f_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &f.sign);
    //    queue.enqueueWriteBuffer(f_buffer, CL_TRUE, 0, sizeof(uint32_t) * f.precision, f.data);
    //    queue.enqueueWriteBuffer(g_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &g.sign);
    //    queue.enqueueWriteBuffer(g_buffer, CL_TRUE, 0, sizeof(uint32_t) * g.precision, g.data);

    //    cl::Kernel kernel(program, "cl_main");
    //    kernel.setArg(0, f_sign_buffer);
    //    kernel.setArg(1, f_buffer);
    //    kernel.setArg(2, g_sign_buffer);
    //    kernel.setArg(3, g_buffer);

    //    cl::Event event;
    //    queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(1), cl::NullRange, nullptr, &event);
    //    event.wait();
    //    queue.enqueueReadBuffer(f_sign_buffer, CL_TRUE, 0, sizeof(int32_t), &f.sign, nullptr, &event);
    //    queue.enqueueReadBuffer(f_buffer, CL_TRUE, 0, sizeof(uint32_t) * f.precision, f.data, nullptr, &event);

    //    std::cout << f.to_fp_string(0x10) << std::endl;
    //}catch(cl::Error err){
    //    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    //}catch(std::exception err){
    //    std::cerr << "ERROR: " << err.what() << std::endl;
    //}

    //return 0;
}
