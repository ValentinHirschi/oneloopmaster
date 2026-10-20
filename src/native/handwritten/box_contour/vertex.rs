//! The contour box's test for the reflection of `zero_external`
//! (`box_complex::zero_external_suspicious`), in a module of its own: added to
//! `box_contour` itself it changed how the crate's codegen units split the
//! kernel's hot functions (not inlined across them) and cost 3.5 % more
//! instructions per D0 call on boxes that never reach it.
use super::*;

/// `box_complex::tfun_vertex_suspicious`: whether the contour that `tfun`
/// selects has a pole whose offset `y + beta` from the vertex `-beta` (the
/// image of the corner (1, 0)) the rounding resolves only to a factor 16 of
/// itself (`max(|beta|, |y|) > 16 |y + beta|`), or, with the pole's own
/// rounding, to 2^12 of itself: `|beta| + dy > 4096 |y + beta|`, with
/// `dy = (Kt |y|² + Lt |y| + Nt)/|2 K y + L|` the root's sensitivity to the
/// rounding of `K y² + L y + N` (`Kt`, `Lt`, `Nt` the magnitudes of the
/// products that form them). The contour is selected as in `tfun`; false
/// for its `t1` branches.
fn tfun_vertex_suspicious(
    [
        a,
        b,
        c,
        gin,
        hin,
        d,
        e,
        f,
        jin,
        dpe,
        _dpf,
        dd_,
        de,
        df,
        djin,
        ddpe,
        _ddpf,
        abc,
    ]: [Z; 18],
) -> bool {
    if !(nz(b) && nz(a)) {
        return false;
    }
    let sj = isign(jin, isign(djin, -1));
    let s = r(-(sj as f64));
    let g = s * gin;
    let h = s * hin;
    let j = s * jin;
    let dj = s * djin;
    let [root1, root2] = roots(b, c, a);
    let [o1, o2] = [root1 - ONE, root2 - ONE];
    let w1 = if negative(norm(o1) - norm(o2)) {
        dv(abc, b * o2)
    } else {
        o1
    };
    let w2 = if negative(norm(o2) - norm(o1)) {
        dv(abc, b * o1)
    } else {
        o2
    };
    let swap = negative(norm(root2) - norm(root1));
    let ((beta1, w1), (beta2, w2)) = if swap {
        ((root2, w2), (root1, w1))
    } else {
        ((root1, w1), (root2, w2))
    };
    let data1 = beta_data(
        beta1, w1, b, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, false,
    );
    let min = |x: f64, y: f64| if negative(x - y) { x } else { y };
    let d1 = min(distance(beta1, &data1.y[0]), distance(beta1, &data1.y[1]));
    // The second contour is needed only where the first one encloses a pole
    // (`d1` not 1), where `tfun` compares their distances.
    let data2;
    let (beta, data) = if d1 - 1.0 != 0.0 {
        data2 = beta_data(
            beta2, w2, b, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, false,
        );
        let d2 = min(distance(beta2, &data2.y[0]), distance(beta2, &data2.y[1]));
        if negative(d2 - d1) {
            (beta1, &data1)
        } else {
            (beta2, &data2)
        }
    } else {
        (beta1, &data1)
    };
    let size = |z: Z| norm(z).sqrt();
    let [a1, b1, c1] = [data.a, data.b, data.c];
    let kt = size(b) * size(a1) + size(h) * size(b1);
    let lt = size(e) * size(a1) + size(h) * size(c1) + size(j) * size(b1);
    let nt = size(f) * size(a1) + size(j) * size(c1);
    let l = e * a1 - h * c1 - j * b1;
    let mb = size(beta);
    data.y.iter().any(|y| {
        let offset = size(y.v + beta);
        let my = size(y.v);
        let larger = if negative(mb - my) { my } else { mb };
        let geometric = negative(16.0 * offset - larger);
        let dy = (kt * my * my + lt * my + nt) / size(r(2.0) * data.k * y.v + l);
        geometric || negative(4096.0 * offset - (mb + dy))
    })
}

/// `box_complex::zero_external_suspicious`: whether one of the two `tfun` of
/// `zero_external_chart` has a selected pole that `tfun_vertex_suspicious`
/// shows to be inaccurate next to the corner (1, 0). Its coefficients are
/// formed as in `zero_external_chart`.
#[cold]
#[inline(never)]
pub(super) fn zero_external_suspicious(p: [Z; 6], m: [Z; 4]) -> bool {
    let dm = mass_directions(m);
    let a = p[2];
    let b = p[1];
    let c = p[5] - p[1] - p[2];
    let h = p[3] - p[4] - p[5] + p[1];
    let j = p[4] - p[0] - p[1];
    let d = m[2] - m[3] - p[2];
    let e = m[1] - m[2] - p[5] + p[2];
    let kk = m[0] - m[1] + p[5] - p[3];
    let f = m[3];
    let cph = p[3] - p[4] - p[2];
    let dpe = m[1] - m[3] - p[5];
    let epk = m[0] - m[2] + p[2] - p[3];
    let dek = m[0] - m[3] - p[3];
    let dpf = m[2] - p[2];
    let bj = p[4] - p[0];
    tfun_vertex_suspicious([
        a,
        b,
        c,
        h,
        j,
        d,
        e,
        f,
        kk,
        dpe,
        dpf,
        dm[2] - dm[3],
        dm[1] - dm[2],
        dm[3],
        dm[0] - dm[1],
        dm[1] - dm[3],
        dm[2],
        p[5],
    ]) || tfun_vertex_suspicious([
        a,
        bj,
        cph,
        h,
        j,
        d,
        epk,
        f,
        kk,
        dek,
        dpf,
        dm[2] - dm[3],
        dm[0] - dm[2],
        dm[3],
        dm[0] - dm[1],
        dm[0] - dm[3],
        dm[2],
        p[3] - p[0],
    ])
}
