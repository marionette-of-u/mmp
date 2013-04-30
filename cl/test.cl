__kernel void test(
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

    div(&h, &f, &g);

    // out f
    *sign_f = h.sign;
    for(base_type i = 0; i < PREC; ++i){
        f_[i] = h.data[i];
    }
}
