#ifndef RK4_CORE_ELEMENT_H
#define RK4_CORE_ELEMENT_H

//----------------------------------------------------------------------
// 化合物中的一个元素项：a=质量数，z=原子序数，s=化学计量数
// 放在 core/ 是因为核素质量表与靶材密度都要用它，
// 避免 core 反过来依赖 sim。
//----------------------------------------------------------------------

struct Element {
    int a;  // 质量数
    int z;  // 原子序数
    double s;  // Relative atomic count; fractional counts represent gas mole mixtures.
};

#endif // RK4_CORE_ELEMENT_H
