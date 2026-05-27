#pragma once

#include "Base.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace Utils {

// 1. 計算二維向量交叉乘積 (Cross Product)，用來判斷旋轉方向
inline float
crossProduct(const Particle &O, const Particle &A, const Particle &B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// 2. Monotone Chain 凸包演算法，時間複雜度 O(K log K)
inline std::vector<Particle> getConvexHull(std::vector<Particle> &pts) {
    int n = pts.size(), k = 0;
    if (n <= 3)
        return pts;
    std::vector<Particle> hull(2 * n);

    // 依照 X 軸排序，若 X 相同則依 Y 軸排序
    std::sort(pts.begin(), pts.end(), [](const Particle &a, const Particle &b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });

    // 建立下凸包 (Lower hull)
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && crossProduct(hull[k - 2], hull[k - 1], pts[i]) <= 0)
            k--;
        hull[k++] = pts[i];
    }

    // 建立上凸包 (Upper hull)
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && crossProduct(hull[k - 2], hull[k - 1], pts[i]) <= 0)
            k--;
        hull[k++] = pts[i];
    }

    hull.resize(k - 1); // 移除重複的起點
    return hull;
}

// 3. 計算點 C 到線段 AB 的最短距離
inline float pointToSegmentDistance(float x, float y, Particle A, Particle B) {
    float ab_x = B.x - A.x;
    float ab_y = B.y - A.y;
    float ac_x = x - A.x;
    float ac_y = y - A.y;

    float ab_sq = ab_x * ab_x + ab_y * ab_y;
    if (ab_sq < 1e-8f) {
        // A 與 B 基本上是同一點
        return std::sqrt(ac_x * ac_x + ac_y * ac_y);
    }

    // 計算投影比例常數 t
    float t = (ac_x * ab_x + ac_y * ab_y) / ab_sq;
    t = std::max(0.0f, std::min(1.0f, t)); // 將投影限制在线段內

    // 計算最近點 Q 的坐標
    float qx = A.x + t * ab_x;
    float qy = A.y + t * ab_y;

    float dx = x - qx;
    float dy = y - qy;
    return std::sqrt(dx * dx + dy * dy);
}

// 4. 主函數：計算圓心 C 到粒子群凸包邊界的最小半徑
inline float getMinRadiusFromHull(
    float x, float y, std::vector<Particle> &particles_in_circle
) {
    if (particles_in_circle.empty())
        return 0.0f;

    // 計算凸包
    std::vector<Particle> hull = getConvexHull(particles_in_circle);

    // 如果凸包頂點太少（退化為點或線段），退化為直接計算到各點的最短距離
    if (hull.size() < 3) {
        float min_dist = 1e9f;
        for (const auto &p : hull) {
            float dx = x - p.x;
            float dy = y - p.y;
            min_dist = std::min(min_dist, std::sqrt(dx * dx + dy * dy));
        }
        return min_dist;
    }

    // 遍歷凸包的每一條邊，計算 C 到線段的最小距離
    float min_radius = 1e9f;
    int m = hull.size();
    for (int i = 0; i < m; i++) {
        Particle A = hull[i];
        Particle B = hull[(i + 1) % m]; // 環形連接相鄰頂點

        float dist = pointToSegmentDistance(x, y, A, B);
        min_radius = std::min(min_radius, dist);
    }

    return min_radius;
}

}; // namespace Utils
