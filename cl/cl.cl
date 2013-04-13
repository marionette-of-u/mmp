// PREC = 1[integral part] + n[fraction part]
#define INTEGER_PART 1
#define FRACTION_PART 4

#define PREC (INTEGER_PART + FRACTION_PART)

typedef uint u_base_type;
typedef int base_type;
typedef ulong u_base2_type;
typedef long base2_type;
#define BASE2_TYPE_SIZE 64
#define BASE2_TYPE_MASK 0x00000000FFFFFFFFUL

typedef u_base_type raw_data_type;

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

// z. := v. * w.
void fp_mul(fixed_point *z, const fixed_point *v, const fixed_point *w){
    if(v->sign == 0 || w->sign == 0){
        z->sign = 0;
        return;
    }
    zero_clear(z->data);
    z->sign = v->sign * w->sign;
    u_base_type buff[PREC * 2];
    for(base_type i = 0; i < PREC * 2; ++i){ buff[i] = 0; }
    for(base_type i = 0; i < PREC; ++i){
        u_base2_type carry2 = 0;
        for(base_type j = 0; j < PREC; ++j){
            base_type k = i + j;
            u_base2_type n = (u_base2_type)(v->data[i]) * (u_base2_type)(w->data[j]) + carry2;
            buff[k] += (u_base_type)(n & BASE2_TYPE_MASK);
            carry2 = (u_base_type)(n >> (BASE2_TYPE_SIZE / 2));
        }
        buff[i + PREC] += (u_base_type)carry2;
    }
    for(base_type i = 0; i < PREC; ++i){
        z->data[i] = buff[PREC + i - 1];
    }
}

// r := z / w
void div_by_word(fixed_point *z, const fixed_point *v, u_base_type w){
    set_zero(z);
    zero_clear(z->data);
    u_base2_type dividend = 0;
    for(base_type i = PREC - 1; i >= 0; --i){
        dividend |= v->data[i];
        z->data[i] = (u_base_type)(dividend / w);
        dividend = (dividend % w) << (BASE2_TYPE_SIZE / 2);
    }
    zero_normalize(z);
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

    fp_mul(&h, &f, &g);

    // out f
    *sign_f = h.sign;
    for(base_type i = 0; i < PREC; ++i){
        f_[i] = h.data[i];
    }
}
