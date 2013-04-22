#include "mmp.hpp"

int main(){
    try{
        // 初期化
        // 整数部精度 : 1 ワード
        // 小数部精度 : 4 ワード
        fixed_point
            f(1, 4, "1", "0"), // 1.0
            g(1, 4, "3", "0"); // 3.0

        // 固定小数点数として表示
        std::cout << f.to_fp_string() << std::endl;
        std::cout << g.to_fp_string() << std::endl;

        // GPU で fixed point を演算するために初期化
        // cl/test.cl を結合して build
        // 整数部精度 : 1 ワード
        // 小数部精度 : 4 ワード
        mmp mmp_manager_test(1, 4, "cl/test.cl");
        {
            // cl/test.cl 内の kernel test 関数で演算
            mmp::kernel_functor functor(mmp_manager_test.create_functor("test"));

            // 引数をセット
            functor.set_arg(0, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, 1, &f.sign);
            functor.set_arg(1, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, f.precision, f.data);
            functor.set_arg(2, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, 1, &g.sign);
            functor.set_arg(3, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, g.precision, g.data);

            // 実行
            // * 結果が返ってくるまで処理は戻らない
            functor.launch(cl::NDRange(1), cl::NullRange);

            // 結果を読む
            functor.read(0, 1, &f.sign);
            functor.read(1, f.precision, f.data);

            // 結果を表示する
            std::cout << f.to_fp_string(0x10) << std::endl;
        }

        // 精度を変更して rebuild する
        // 整数部精度 : 1 ワード
        // 小数部精度 : 8 ワード
        mmp_manager_test.rebuild(1, 8, "cl/test.cl");
        f.set_prec(1, 8);
        f.set("1", "0");
        g.set_prec(1, 8);
        g.set("3", "0");
        {
            mmp::kernel_functor functor(mmp_manager_test.create_functor("test"));
            functor.set_arg(0, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, 1, &f.sign);
            functor.set_arg(1, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, f.precision, f.data);
            functor.set_arg(2, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, 1, &g.sign);
            functor.set_arg(3, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, g.precision, g.data);
            functor.launch(cl::NDRange(1), cl::NullRange);
            functor.read(0, 1, &f.sign);
            functor.read(1, f.precision, f.data);
            std::cout << f.to_fp_string(0x10) << std::endl;
        }
    }catch(cl::Error err){
        std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    }catch(std::exception err){
        std::cerr << "ERROR: " << err.what() << std::endl;
    }

    return 0;
}
