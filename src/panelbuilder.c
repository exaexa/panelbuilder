
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

void
pw_gd(const int* np,
      const int* dp,
      const int* kp,
      const int* itersp,
      const float* alphap,
      const float* tolp,
      const float* s_dk,
      const float* spw_dk,
      const float* snw_dk,
      const float* nw_k,
      const float* y_dn,
      float* x_kn,
      float* r_dn,
      float* g_k)
{
  const size_t n = *np, d = *dp, k = *kp, iters = *itersp;
  const float alpha = *alphap, tol = *tolp;
  size_t ni, di, ki, ii;

  /* x_kn are pre-initialized to starting points */
  for (ni = 0; ni < n; ++ni) {
    float* restrict x_k = x_kn + ni * k;
    float* restrict r_d = r_dn + ni * d;
    const float* restrict y_d = y_dn + ni * d;

    for (ii = 0; ii < iters; ++ii) {
      /* compute the residuals */
      for (di = 0; di < d; ++di)
        r_d[di] = -y_d[di];
      for (ki = 0; ki < k; ++ki)
        for (di = 0; di < d; ++di)
          r_d[di] += x_k[ki] * s_dk[di + d * ki];

      float acc = 0;
      for (ki = 0; ki < k; ++ki) {
        /* guess the direction */
        float gki = (x_k[ki] > 0 ? 0 : nw_k[ki] * x_k[ki]);
        for (di = 0; di < d; ++di) {
          const float w =
            r_d[di] > 0 ? spw_dk[di + d * ki] : snw_dk[di + d * ki];
          gki += r_d[di] * s_dk[di + d * ki] * w;
        }

        /* apply the gradient */
        x_k[ki] -= alpha * gki;
        acc += gki * gki;
      }

      /* if we didn't move too much, terminate */
      if (acc < tol)
        break;
    }
  }
}

static const R_CMethodDef cMethods[] = { { "pw_gd", (DL_FUNC)&pw_gd, 14 },
                                         { NULL, NULL, 0 } };

void
R_init_panelbuilder(DllInfo* info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
