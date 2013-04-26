#ifndef MMP_HPP_
#define MMP_HPP_

#include <algorithm>
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <list>
#include <limits>
#include <sstream>
#include <string>
#include <exception>
#include <cmath>
#include <cstdint>
#include <cstring>
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#define BASE2_TYPE_SIZE 64
#define BASE2_TYPE_MASK 0x00000000FFFFFFFFUL

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

template<class Signature = void>
class x_fixed_point{
private:
    std::size_t m_integral_part, m_fraction_part, m_precision;

public:
    const std::size_t &integral_part, &fraction_part, &precision;
    int32_t sign;
    uint32_t *data;

public:
    x_fixed_point(std::size_t integral_part_, std::size_t fraction_part_, int n) :
        m_integral_part(integral_part_),
        m_fraction_part(fraction_part_),
        m_precision(integral_part_ + fraction_part_),
        integral_part(m_integral_part),
        fraction_part(m_fraction_part),
        precision(m_precision),
        sign(0),
        data(new uint32_t[precision])
    { set(n); }

    x_fixed_point(
        std::size_t integral_part_,
        std::size_t fraction_part_,
        const std::string &str_i,
        const std::string &str_f,
        uint32_t radix = 10
    ) :
        m_integral_part(integral_part_),
        m_fraction_part(fraction_part_),
        m_precision(integral_part_ + fraction_part_),
        integral_part(m_integral_part),
        fraction_part(m_fraction_part),
        precision(m_precision),
        sign(0),
        data(new uint32_t[precision])
    { set(str_i, str_f, radix); }

    x_fixed_point(const x_fixed_point &other) :
        m_integral_part(0),
        m_fraction_part(0),
        m_precision(0),
        integral_part(m_integral_part),
        fraction_part(m_fraction_part),
        precision(m_precision),
        sign(0),
        data(new uint32_t[other.precision])
    { set(other); }

    x_fixed_point(x_fixed_point &&other) :
        m_integral_part(other.m_integral_part),
        m_fraction_part(other.m_fraction_part),
        m_precision(other.m_precision),
        integral_part(m_integral_part),
        fraction_part(m_fraction_part),
        precision(m_precision),
        sign(other.sign),
        data(other.data)
    { other.data = nullptr; }

    ~x_fixed_point(){
        delete[] data;
    }

    std::string to_string(uint32_t radix = 10) const{
        if(sign == 0){ return "0"; }
        std::string result;
        x_fixed_point temp(integral_part, fraction_part, 0);
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
            x_fixed_point integer(integral_part, 0, 0);
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
        x_fixed_point temp(integral_part, fraction_part, 0);
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

    void set(int n){
        sign = n > 0 ? 1 : n < 0 ? -1 : 0;
        for(size_t i = 0; i < precision; ++i){ data[i] = 0; }
        data[fraction_part] = static_cast<uint32_t>(n < 0 ? -n : n);
    }

    void set(int n, const std::string &str_i, const std::string &str_f, uint32_t radix = 10){
        set(str_i, str_f, radix);
        sign = n > 0 ? 1 : n < 0 ? -1 : 0;
        if(sign == 0){ primitive_set_zero(); }
    }

    void set(const std::string &str_i, const std::string &str_f, uint32_t radix = 10){
        sign = 1;
        for(std::size_t i = 0; i < precision; ++i){ data[i] = 0; }
        {
            std::size_t n = static_cast<std::size_t>(std::ceil((std::log(static_cast<double>(radix)) / std::log(2.0)) * str_i.size() / (BASE2_TYPE_SIZE / 2)));
            if(n > integral_part){ throw(fixed_point_exception("input value is too large.")); }
            x_fixed_point digit(integral_part, fraction_part, 1), z(integral_part, fraction_part, 1);
            for(size_t i = 0; i < precision; ++i){ z.data[i] = 0; }
            char temp[] = { 0, '\0' };
            for(std::size_t i = 0; i < str_i.size(); ++i){
                z.primitive_mul_by_single(radix);
                temp[0] = str_i[i];
                digit.data[fraction_part] = std::strtol(temp, nullptr, radix);
                z.primitive_add_n(digit);
            }
            primitive_add_n(z);
        }
        {
            x_fixed_point digit(integral_part, fraction_part, 1), z(integral_part, fraction_part, 1);
            for(size_t i = 0; i < precision; ++i){ z.data[i] = 0; }
            digit.data[fraction_part] = 0;
            char temp[] = { 0, '\0' };
            for(int i = static_cast<int>(str_f.size() - 1); i >= 0; --i){
                temp[0] = str_f[i];
                digit.data[fraction_part] = std::strtol(temp, nullptr, radix);
                z.primitive_add_n(digit);
                z.primitive_div_by_single(radix);
            }
            primitive_add_n(z);
        }
    }

    void set(const x_fixed_point &w){
        primitive_set_fixed_point(w);
    }

    void set_prec(std::size_t i_part, std::size_t f_part){
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
        m_fraction_part = f_part;
        m_integral_part = i_part;
        m_precision = fraction_part + integral_part;
        delete[] data;
        data = new_data;
    }

    void mul(uint32_t w){
        primitive_mul_by_single(w);
    }

    void mul(const x_fixed_point &v, const x_fixed_point &w){
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

    void div(uint32_t w){
        primitive_div_by_single(w);
    }

    void div(const x_fixed_point &u, const x_fixed_point &v){
        if(u.sign == 0){
            primitive_set_zero();
            return;
        }
        sign = u.sign * v.sign;
        const int32_t m = u.primitive_degree(), n = v.primitive_degree();
        if(n > m){
            sign = 0;
            primitive_set_zero();
            return;
        }
        const uint64_t b = static_cast<uint64_t>(BASE2_TYPE_MASK) + 1;
        uint32_t *un, *vn;
        uint64_t qhat, rhat, p;
        int64_t i, j, t, k, z;
        for(i = 0; i < static_cast<int64_t>(precision); ++i){ data[i] = 0; }
        if(n == 1){
            k = 0;
            for(j = m - 1; j >= 0; --j){
                z = k * b + u.data[j];
                data[j] = static_cast<uint32_t>(z / v.data[0]);
                k = z - data[j] * v.data[0];
            }
            return;
        }
        const uint32_t s = primitive_nlz(v.data[n - 1]);
        std::unique_ptr<uint32_t[]> vn_scoped_guard(new uint32_t[precision * 2]);
        vn = vn_scoped_guard.get();
        for(i = n - 1; i > 0; --i){
            vn[i] = static_cast<uint32_t>(primitive_s_lshift(v.data[i], s) | primitive_s_rshift(v.data[i - 1], BASE2_TYPE_SIZE / 2 - s));
        }
        vn[0] = static_cast<uint32_t>(primitive_s_lshift(v.data[0], s));
        std::unique_ptr<uint32_t[]> un_scoped_guard(new uint32_t[(precision + 1) * 2]);
        un = un_scoped_guard.get();
        un[m + fraction_part] = static_cast<uint32_t>(primitive_s_rshift(u.data[m - 1], BASE2_TYPE_SIZE / 2 - s));
        for(i = 0; i < static_cast<int64_t>(precision); ++i){ un[i] = 0; }
        for(i = m - 1; i > 0; --i){
            un[i + fraction_part] = static_cast<uint32_t>(primitive_s_lshift(u.data[i], s) | primitive_s_rshift(u.data[i - 1], BASE2_TYPE_SIZE / 2 - s));
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

    void add(const x_fixed_point &w){
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
                x_fixed_point temp(integral_part, fraction_part, 1);
                temp.primitive_set_fixed_point(*this);
                primitive_set_fixed_point(w);
                primitive_sub_n(temp);
            }
            primitive_zero_normalize();
        }else{
            primitive_set_zero();
        }
    }

    void sub(const x_fixed_point &w){
        if(w.sign == 0){ return; }
        if(sign == 0){
            primitive_set_fixed_point(w);
            sign = -sign;
            return;
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
                x_fixed_point temp(integral_part, fraction_part, 1);
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

    static int compare(const x_fixed_point &x, const x_fixed_point &y){
        int n = primitive_compare(x, y);
        if(n == 0){
            return x.sign > y.sign ? 1 : x.sign < y.sign ? -1 : 0;
        }else{
            return n;
        }
    }

private:
    void primitive_mul_inv(){
        if((data[0] & 1) == 0){ data[0] += 1; }
        x_fixed_point
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

    void primitive_mul_with_sign(const x_fixed_point &v, const x_fixed_point &w){
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

    void primitive_add_n(const x_fixed_point &w){
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

    void primitive_sub_n(const x_fixed_point &w){
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
        x_fixed_point result(integral_part, fraction_part, 1);
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

    void primitive_set_fixed_point(const x_fixed_point &w){
        if(this == &w){ return; }
        m_integral_part = w.integral_part;
        m_fraction_part = w.fraction_part;
        m_precision = w.precision;
        sign = w.sign;
        for(std::size_t i = 0; i < precision; ++i){ data[i] = w.data[i]; }
    }

    static int primitive_compare(const x_fixed_point &z, const x_fixed_point &w){
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

typedef x_fixed_point<> fixed_point;

template<class Signature = void>
class x_mmp{
    friend class kernel_functor;

private:
    struct kernel_functor_parameter{
        cl::Program *program;
        std::string fn_name;
        cl::Context *context;
        cl::CommandQueue *queue;
    };

public:
    class kernel_functor{
    public:
        struct mem_flag{
            enum{
                rw = CL_MEM_READ_WRITE,
                w = CL_MEM_WRITE_ONLY,
                r = CL_MEM_READ_ONLY,
                use_host_ptr = CL_MEM_USE_HOST_PTR,
                alloc_host_ptr = CL_MEM_ALLOC_HOST_PTR,
                copy_host_ptr = CL_MEM_COPY_HOST_PTR
            };
        };

        kernel_functor(kernel_functor_parameter param) :
            buffer_list(),
            kernel_instance(*param.program, param.fn_name.c_str()),
            event(nullptr),
            context_ptr(param.context),
            queue_ptr(param.queue)
        {}

        template<class T>
        void set_arg(cl_uint idx, const T &value){
            kernel_instance.setArg(idx, value);
        }

        template<class T>
        void set_arg(cl_uint idx, cl_mem_flags flags, std::size_t data_size, T *data){
            cl::Buffer *buffer = new cl::Buffer(*context_ptr, flags, sizeof(T) * data_size, static_cast<void*>(data));
            buffer_list[idx].reset(buffer);
            queue_ptr->enqueueWriteBuffer(*buffer, CL_TRUE, 0, data_size, data);
            kernel_instance.setArg(idx, *buffer);
        }

        void launch(const cl::NDRange &global, const cl::NDRange &local){
            event.reset(new cl::Event);
            queue_ptr->enqueueNDRangeKernel(kernel_instance, cl::NullRange, global, local, nullptr, event.get());
            event->wait();
        }

        template<class T>
        void read(cl_uint idx, std::size_t data_size, T *data){
            queue_ptr->enqueueReadBuffer(*(buffer_list[idx]), CL_TRUE, 0, sizeof(T) * data_size, data, nullptr);
        }

    private:
        std::map<cl_uint, std::unique_ptr<cl::Buffer>> buffer_list;
        cl::Kernel kernel_instance;
        std::unique_ptr<cl::Event> event;
        cl::Context *context_ptr;
        cl::CommandQueue *queue_ptr;

        kernel_functor(){}
    };

    x_mmp(std::size_t i_part, std::size_t f_part, std::string mmp_kernel_filepath, std::string mmp_filepath_ = "cl/mmp.cl") :
        mmp_filepath(mmp_filepath_),
        integral_part(i_part),
        fraction_part(f_part),
        platform(),
        context(nullptr),
        devices(nullptr),
        queue(nullptr),
        sources(nullptr),
        program(nullptr)
    {
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if(platforms.size() == 0){
            throw(std::exception("platform size 0."));
        }
        platform = platforms[0];
        cl_context_properties properties[] = {
            CL_CONTEXT_PLATFORM,
            reinterpret_cast<cl_context_properties>(platform()),
            0
        };
        context.reset(new cl::Context(CL_DEVICE_TYPE_GPU, properties));
        devices.reset(new std::vector<cl::Device>);
        *devices = context->getInfo<CL_CONTEXT_DEVICES>();
        queue.reset(new cl::CommandQueue(*context, (*devices)[0]));
        rebuild(integral_part, fraction_part, mmp_kernel_filepath);
    }

    std::size_t get_integral_part() const{
        return integral_part;
    }

    std::size_t get_fraction_part() const{
        return fraction_part;
    }

    kernel_functor_parameter create_functor(const std::string &fn_name){
        kernel_functor_parameter param;
        param.program = program.get();
        param.fn_name = fn_name;
        param.context = context.get();
        param.queue = queue.get();
        return param;
    }

    void rebuild(std::size_t i_part, std::size_t f_part, const std::string &mmp_kernel_filepath){
        raw_sources.clear();
        sources.reset(new cl::Program::Sources);
        std::stringstream ss;
        ss << "#define INTEGER_PART " << i_part << "\n";
        ss << "#define FRACTION_PART " << f_part << "\n";
        std::string i_f_part = ss.str();
        raw_sources.insert(raw_sources.begin(), i_f_part.begin(), i_f_part.end());
        push_back_source(mmp_filepath);
        push_back_source(mmp_kernel_filepath);
        sources->push_back(std::make_pair(&raw_sources[0], raw_sources.size()));
        program.reset(new cl::Program(*context, *sources));
        program->build(*devices, "-cl-strict-aliasing");
    }

private:
    void push_back_source(const std::string &filepath){
        std::vector<char> cl_source;
        std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::ate);
        if(ifile.fail()){ throw(std::exception("can not open cl file.")); }
        cl_source.resize(static_cast<std::size_t>(ifile.tellg()));
        ifile.seekg(0, ifile.beg);
        ifile.read(&cl_source[0], cl_source.size());
        raw_sources.insert(raw_sources.end(), cl_source.begin(), cl_source.end());
    }

    const std::string mmp_filepath;
    std::size_t integral_part, fraction_part;
    cl::Platform platform;
    std::vector<char> raw_sources;
    std::unique_ptr<cl::Context> context;
    std::unique_ptr<std::vector<cl::Device>> devices;
    std::unique_ptr<cl::CommandQueue> queue;
    std::unique_ptr<cl::Program::Sources> sources;
    std::unique_ptr<cl::Program> program;
};

typedef x_mmp<> mmp;

#endif // MMP_HPP_
