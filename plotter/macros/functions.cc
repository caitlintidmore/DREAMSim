#include "TMath.h"

double CosTheta(double x0, double y0, double z0, double px, double py, double pz)
{
    double c0 = 30.0 / 1.49; // cm / ns
    double p = TMath::Sqrt(px * px + py * py + pz * pz);
    px = px / p * c0;
    py = py / p * c0;
    pz = pz / p * c0;
    p = TMath::Sqrt(px * px + py * py + pz * pz);

    double r0 = 0.039;
    double a = px * px + py * py;
    double b = 2 * (x0 * px + y0 * py);
    double c = x0 * x0 + y0 * y0 - r0 * r0;

    if (b * b - 4 * a * c < 0)
    {
        return -2.0;
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
        return -1.0;
    }

    double x1 = x0 + t * px;
    double y1 = y0 + t * py;
    double z1 = z0 + t * pz;

    return (px * x1 + py * y1) / (p * r0);
}

double SinTheta(double x0, double y0, double z0, double px, double py, double pz)
{
    double costheta = CosTheta(x0, y0, z0, px, py, pz);
    if (costheta < 0)
    {
        return costheta;
    }
    return TMath::Sqrt(1 - costheta * costheta);
}