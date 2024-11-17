#include "TMath.h"

#include "ROOT/RVec.hxx"
using namespace ROOT::VecOps;

RVec<double> CosTheta(const RVec<double> &x0, const RVec<double> &y0, const RVec<double> &z0, const RVec<double> &px, const RVec<double> &py, const RVec<double> &pz)
{
    RVec<double> result(x0.size());

    double c0 = 30.0 / 1.49; // cm / ns
    double r0 = 0.039;       // cm

    for (size_t i = 0; i < x0.size(); ++i)
    {
        double p = TMath::Sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
        double px_norm = px[i] / p * c0;
        double py_norm = py[i] / p * c0;
        double pz_norm = pz[i] / p * c0;
        p = TMath::Sqrt(px_norm * px_norm + py_norm * py_norm + pz_norm * pz_norm);

        double a = px_norm * px_norm + py_norm * py_norm;
        double b = 2 * (x0[i] * px_norm + y0[i] * py_norm);
        double c = x0[i] * x0[i] + y0[i] * y0[i] - r0 * r0;

        if (b * b - 4 * a * c < 0)
        {
            result[i] = -2.0;
            continue;
        }

        double t1 = (-b + TMath::Sqrt(b * b - 4 * a * c)) / (2 * a); // ns
        double t2 = (-b - TMath::Sqrt(b * b - 4 * a * c)) / (2 * a); // ns

        double tmin = TMath::Min(t1, t2);
        double tmax = TMath::Max(t1, t2);
        double t = tmin;
        if (tmin < 0)
        {
            t = tmax;
        }
        if (t < 0)
        {
            result[i] = -1.0;
            continue;
        }

        double x1 = x0[i] + t * px_norm;
        double y1 = y0[i] + t * py_norm;
        double z1 = z0[i] + t * pz_norm;

        result[i] = (px_norm * x1 + py_norm * y1) / (p * r0);
    }

    return result;
}

RVec<double> SinTheta(const RVec<double> &x0, const RVec<double> &y0, const RVec<double> &z0, const RVec<double> &px, const RVec<double> &py, const RVec<double> &pz)
{
    RVec<double> costheta = CosTheta(x0, y0, z0, px, py, pz);
    RVec<double> result(costheta.size());

    for (size_t i = 0; i < costheta.size(); ++i)
    {
        if (costheta[i] < 0)
        {
            result[i] = costheta[i];
        }
        else
        {
            result[i] = TMath::Sqrt(1 - costheta[i] * costheta[i]);
        }
    }

    return result;
}