// PREC = 1[integral part] + n[fraction part]
#define PREC 5

typedef uint u_base_type;
typedef int base_type;
typedef ulong u_base2_type;
typedef long base2_type;
#define BASE2_TYPE_SIZE 64
#define BASE2_TYPE_MASK 0x00000000FFFFFFFFUL

typedef u_base_type raw_data_type;

//base_type msb(u_base_type x){
//    if(x == 0){ return -1; }
//    base_type m = 0;
//    if(x & 0xFFFF0000){ m |= 16; x >>= 16; }
//    if(x & 0xFF00){ m |= 8; x >>= 8; }
//    if(x & 0xF0){ m |= 4; x >>= 4; }
//    if(x & 0xC){ m |= 2; x >>= 2; }
//    if(x & 0x2){ m |= 1; }
//    return m;
//}
//
//base_type msbs(const raw_data_type *x){
//    for(base_type i = 0; i < PREC; ++i){
//        base_type m = msb(x[i]);
//        if(m >= 0){ return m + ((PREC - i - 1) << 5); }
//    }
//    return -1;
//}

// compare
base_type compare(const raw_data_type *z, const raw_data_type *w){
    for(base_type i = PREC - 1; i >= 0; --i){
        if(z[i] > w[i]){ return 1; }
        if(z[i] < w[i]){ return -1; }
    }
    return 0;
}

// zero-clear z.
void zero_clear(raw_data_type *z){
    for(base_type i = 0; i < PREC; ++i){
        z[i] = 0;
    }
}

// copy src to dst of n words.
void copy_array(raw_data_type *dst, const raw_data_type *src){
    for(base_type i = 0; i < PREC; ++i){
        dst[i] = src[i];
    }
}

// add n words with carry.
// return value is overflow value.
u_base_type add_n(raw_data_type *z, const raw_data_type *w, u_base_type carry){
    u_base2_type carry2 = carry;
    for(base_type i = 0; i < PREC; ++i){
        carry2 += (u_base2_type)z[i] + (u_base2_type)w[i];
        z[i] = (u_base_type)(carry2 & BASE2_TYPE_MASK);
        carry2 >>= (BASE2_TYPE_SIZE / 2);
    }
    return (u_base_type)carry2;
}

// sub n words with borrow.
// return value is underflow value.
u_base_type sub_n(raw_data_type *z, const raw_data_type *w, u_base_type borrow){
    u_base2_type borrow2 = borrow;
    for(base_type i = 0; i < PREC; ++i){
        u_base2_type x = (u_base2_type)z[i] - (u_base2_type)w[i] - borrow2;
        z[i] = (u_base_type)(x & BASE2_TYPE_MASK);
        borrow2 = x > BASE2_TYPE_MASK ? 1 : 0;
    }
    return (u_base_type)borrow2;
}

// mul multiprec by singleprec with carry.
// return value is overflow value.
u_base_type mul_by_single(raw_data_type *z, u_base_type w, u_base_type carry){
    u_base2_type carry2 = carry;
    for(base_type i = 0; i < PREC; ++i){
        u_base2_type x = (u_base2_type)z[i] * (u_base2_type)w + (u_base2_type)carry2;
        z[i] = (u_base_type)(x & BASE2_TYPE_MASK);
        carry2 = x >> (u_base2_type)(BASE2_TYPE_SIZE / 2);
    }
    return (u_base_type)carry2;
}

// div multiprec by singleprec with remainder.
// return value is remainder value.
u_base_type div_by_single(raw_data_type *z, u_base_type w, u_base_type rem){
    u_base2_type r = rem;
    for(base_type i = PREC - 1; i >= 0; --i){
        u_base2_type x = (r << (u_base2_type)(BASE2_TYPE_SIZE / 2)) | (u_base2_type)z[i];
        z[i] = (base_type)(x / (u_base2_type)w);
        r = x % (u_base2_type)w;
    }
    return (u_base_type)r;
}

u_base_type shift_left(raw_data_type *z, u_base_type bits, u_base_type k){
    u_base2_type mask = ((u_base2_type)1 << (u_base2_type)k) - 1;
    u_base2_type carry = (u_base2_type)bits & mask;
    for(base_type i = 0; i < PREC; ++i){
        u_base2_type x = ((u_base2_type)z[i] << (u_base2_type)k) | carry;
        z[i] = (u_base_type)(x & BASE2_TYPE_MASK);
        carry = x >> (u_base2_type)(BASE2_TYPE_SIZE / 2);
    }
    return (u_base_type)carry;
}

u_base_type shift_right(raw_data_type *z, u_base_type bits, u_base_type k){
    u_base2_type mask = ((u_base2_type)1 << (u_base2_type)k) - 1;
    u_base2_type carry = (u_base2_type)bits & mask;
    for(base_type i = PREC - 1; i >= 0; --i){
        u_base2_type x = (carry << (u_base2_type)(BASE2_TYPE_SIZE / 2)) | (u_base2_type)z[i];
        z[i] = (u_base_type)(x >> (u_base2_type)k);
        carry = x & mask;
    }
    return (u_base_type)carry;
}

bool check_zero(const raw_data_type *z){
    for(base_type i = 0; i < PREC; ++i){
        if(z[i] != 0){ return false; }
    }
    return true;
}

typedef struct{
    base_type sign;
    raw_data_type data[PREC];
} fixed_point;

void copy_from_pre_fixed_point(fixed_point *z, const u_base_type *w){
    if(check_zero(w)){
        zero_clear(z);
        return;
    }
    z->sign = 1;
    for(base_type i = 0; i < PREC; ++i){
        z->data[i] = w[i];
    }
}

base_type copy_to_pre_fixed_point(u_base_type *z, const fixed_point *w){
    for(base_type i = 0; i < PREC; ++i){
        z[i] = w->data[i];
    }
    return w->sign;
}

void to_abs(fixed_point *z){
    if(z->sign < 0){
        z->sign = -z->sign;
    }
}

// z := 0
void set_zero(fixed_point *z){
    z->sign = 0;
    //zero_clear(z->data);
}

// z := x
void set_base(fixed_point *z, base_type x){
    set_zero(z);
    if(x == 0){
        z->sign = 0;
        return;
    }
    if(x > 0){
        z->sign = 1;
        z->data[0] = (raw_data_type)x;
    }else{
        z->sign = -1;
        (raw_data_type)-x;
    }
}

// z := x
void set_fixed_point(fixed_point *z, const fixed_point *x){
    if(z == x){ return; }
    z->sign = x->sign;
    copy_array(z->data, x->data);
}

void zero_normalize(fixed_point *z){
    if(z->sign == 0){ return; }
    if(check_zero(z->data)){ z->sign = 0; }
}

void negate(fixed_point *z){
    z->sign = -z->sign;
}

// z += w
void add(fixed_point *z, const fixed_point *w){
    if(w->sign == 0){ return; }
    if(z->sign == 0){
        set_fixed_point(z, w);
        return;
    }
    if(z->sign == w->sign){
        add_n(z->data, w->data, 0);
        return;
    }
    base_type u = compare(z->data, w->data);
    if(u != 0){
        if(u > 0){
            sub_n(z->data, w->data, 0);
        }else{
            raw_data_type temp[PREC];
            copy_array(temp, z->data);
            set_fixed_point(z, w);
            sub_n(z->data, temp, 0);
        }
        zero_normalize(z);
    }else{
        set_zero(z);
    }
}

// z -= w
void sub(fixed_point *z, const fixed_point *w){
    if(w->sign == 0){ return; }
    if(z->sign == 0){
        set_fixed_point(z, w);
        z->sign = -z->sign;
        return;
    }
    if(z->sign != w->sign){
        add_n(z->data, w->data, 0);
        return;
    }
    base_type u = compare(z->data, w->data);
    if(u != 0){
        if(u > 0){
            sub_n(z->data, w->data, 0);
        }else{
            raw_data_type temp[PREC];
            copy_array(temp, z->data);
            set_fixed_point(z, w);
            sub_n(z->data, temp, 0);
        }
        zero_normalize(z);
    }else{
        set_zero(z);
    }
}

// z *= w
void mul_by_word(fixed_point *z, base_type w){
    if(w == 0){
        set_zero(z);
        return;
    }
    if(z->sign == 0){ return; }
    if(w < 0){
        negate(z);
        w = -w;
    }
    if(w == 1){ return; }
    mul_by_single(z->data, (u_base_type)w, 0);
}

// z *= 2^w
u_base_type mul_pow2(fixed_point *z, base_type w){
    if(w == 0){ return 0; }
    if(z->sign == 0){ return 0; }
    if(w > 0){
        u_base_type carry = shift_left(z->data, 0, (u_base_type)w);
        return carry;
    }else{
        u_base_type borrow = shift_right(z->data, 0, (u_base_type)(-w));
        zero_normalize(z);
        return borrow;
    }
}

// z := v * w
void mul(fixed_point *z, const fixed_point *v, const fixed_point *w){
    zero_clear(z->data);
    if(v->sign == 0 || w->sign == 0){
        z->sign = 0;
        return;
    }
    z->sign = v->sign * w->sign;
    u_base2_type carry2 = 0;
    for(base_type i = 0; i < PREC; ++i){
        for(base_type j = 0; j < PREC; ++j){
            base_type k = i + j;
            if(k >= PREC){ continue; }
            u_base2_type n = (u_base2_type)(v->data[i]) * (u_base2_type)(w->data[j]) + carry2;
            z->data[k] += (u_base_type)(n & BASE2_TYPE_MASK);
            carry2 = (u_base_type)(n >> (BASE2_TYPE_SIZE / 2));
        }
    }
}

// z mod w
u_base_type mod_by_word(const fixed_point *z, u_base_type w){
    u_base2_type p = 1;
    u_base2_type r = 0;
    for(base_type i = 0; i < PREC; ++i){
        for(base_type bit = 0; bit < BASE2_TYPE_SIZE / 2; ++bit){
            if((z->data[i] & (1 << bit)) != 0){
                r += p;
                if(r >= w){
                    r -= w;
                }
            }
            p <<= 1;
            if(p >= w){
                p -= w;
            }
        }
    }
    return (u_base_type)r;
}

// r := z mod w
void mod(fixed_point *r, const fixed_point *z, const fixed_point *w){
    fixed_point p_, *p = &p_;
    set_zero(p);
    set_zero(r);
    p->data[0] = 1;
    for(base_type i = 0; i < PREC; ++i){
        for(base_type bit = 0; bit < PREC; ++bit){
            if((z->data[i] & (1 << bit)) != 0){
                add(r, p);
                if(compare(r, w) >= 0){
                    sub(r, w);
                }
            }
            u_base_type carry = mul_pow2(p, 1);
            if(carry > 0 || compare(p, w) >= 0){
                fixed_point temp = *w;
                sub(&temp, p);
                *p = temp;
            }
        }
    }
}

// r := z / w
void div_by_word(fixed_point *r, const fixed_point *z, u_base_type w){
    set_base(r, 1);
    u_base2_type dividend = 0;
    for(base_type i = PREC - 1; i >= 0; --i){
        dividend |= z->data[i];
        r->data[i] = (u_base_type)(dividend / w);
        dividend = (dividend % w) << (BASE2_TYPE_SIZE / 2);
    }
}

void str_to_fixed_point(fixed_point *z, const char *str, u_base_type n){
    fixed_point p_, *p = &p_;
    set_zero(z);
    set_base(p, 1);
    for(base_type i = 0; str[i] != '\0'; ++i){
        mul_by_word(z, n);
        p->data[0] = str[i] - '0';
        add(z, p);
    }
}

void digit_to_fixed_point(fixed_point *z, const char *str){
    str_to_fixed_point(z, str, 10);
}

void hex_to_fixed_point(fixed_point *z, const char *str){
    str_to_fixed_point(z, str, 0x10);
}

void fixed_point_to_str(char *buff, const fixed_point *z_prime, u_base_type n){
    fixed_point z_, *z = &z_;
    fixed_point p_, *p = &p_;
    z_ = *z_prime;
    set_base(p, n);
    if(check_zero(z)){
        buff[0] = '0';
        buff[1] = '\0';
        return;
    }else{
        base_type i;
        for(base_type i = 0; check_zero(z); ++i){
            buff[i] = (char)(mod_by_word(z, n) + '0');
            div_by_word(z, z, n);
        }
        buff[i] = '\0';
    }
}

void fixed_point_to_digit(char *buff, const fixed_point *z){
    fixed_point_to_str(buff, z, 10);
}

void fixed_point_to_hex(char *buff, const fixed_point *z){
    fixed_point_to_str(buff, z, 0x10);
}

__kernel void cl_main(
    __global base_type *sign_f,
    __global u_base_type *f_,
    __constant base_type *sign_g,
    __constant u_base_type *g_
){
    fixed_point f, g, h;

    // init f
    f.sign = *sign_f;
    for(base_type i = 0; i < PREC; ++i){
        f.data[i] = f_[i];
    }
    
    // init g
    g.sign = *sign_g;
    for(base_type i = 0; i < PREC; ++i){
        g.data[i] = g_[i];
    }

    mul(&h, &f, &g);

    // out f
    *sign_f = h.sign;
    for(base_type i = 0; i < PREC; ++i){
        f_[i] = h.data[i];
    }
}
