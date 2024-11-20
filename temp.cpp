static mat3 cross_product_squared_transpose(vec3 const& v)
    {
        auto const a = v[0];
        auto const b = v[1];
        auto const c = v[2];
        auto const a2 = a * a;
        auto const b2 = b * b;
        auto const c2 = c * c;

        mat3 M;

        math::get(M, 0, 0) = b2 + c2;
        math::get(M, 1, 1) = a2 + c2;
        math::get(M, 2, 2) = a2 + b2;

        math::get(M, 1, 0) = -a * b;
        math::get(M, 2, 0) = -a * c;
        math::get(M, 2, 1) = -b * c;

        math::get(M, 0, 1) = -a * b;
        math::get(M, 0, 2) = -a * c;
        math::get(M, 1, 2) = -b * c;

        return M;
    }

[[nodiscard]] static quadric probabilistic_triangle_quadric(pos3 const& mean_p, //
                                                                pos3 const& mean_q,
                                                                pos3 const& mean_r,
                                                                scalar_t stddev)
    {
        auto const sigma = stddev * stddev;
        auto const p = math::to_vec(mean_p);
        auto const q = math::to_vec(mean_q);
        auto const r = math::to_vec(mean_r);

        auto const pxq = math::cross(p, q);
        auto const qxr = math::cross(q, r);
        auto const rxp = math::cross(r, p);

        auto const det_pqr = math::dot(pxq, r);

        auto const cross_pqr = pxq + qxr + rxp;

        auto const pmq = p - q;
        auto const qmr = q - r;
        auto const rmp = r - p;

        mat3 A = math::self_outer_product(cross_pqr) +         //
                 (math::cross_product_squared_transpose(pmq) + //
                  math::cross_product_squared_transpose(qmr) + //
                  math::cross_product_squared_transpose(rmp))
                     * sigma;

        auto ss = sigma * sigma;
        auto ss6 = 6 * ss;
        auto ss2 = 2 * ss;
        math::get(A, 0, 0) += ss6;
        math::get(A, 1, 1) += ss6;
        math::get(A, 2, 2) += ss6;

        vec3 b = cross_pqr * det_pqr;

        b = b - (math::cross(pmq, pxq) + math::cross(qmr, qxr) + math::cross(rmp, rxp)) * sigma;

        b = b + (p + q + r) * ss2;

        scalar_t c = det_pqr * det_pqr;

        c += sigma * (math::dot(pxq, pxq) + math::dot(qxr, qxr) + math::dot(rxp, rxp)); // 3x (a x b)^T M_c (a x b)

        c += ss2 * (math::dot(p, p) + math::dot(q, q) + math::dot(r, r)); // 3x a^T Ci[S_b, S_c] a

        c += ss6 * sigma; // Tr[S_r Ci[S_p, S_q]]

        return quadric::from_coefficients(A, b, c);
    }
