#pragma warning disable RECS0018 // Comparison of floating point numbers with equality operator

using System;

namespace Qef
{
    public struct Vec3
    {
        public float x;
        public float y;
        public float z;

        public Vec3(float x, float y, float z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public Vec3(float f)
            : this(f, f, f)
        { }

        public float this[int i]
        {
            get
            {
                switch (i)
                {
                    case 0:
                        return x;

                    case 1:
                        return y;

                    case 2:
                        return z;

                    default:
                        throw new IndexOutOfRangeException();
                }
            }

            set
            {
                switch (i)
                {
                    case 0:
                        x = value;
                        break;

                    case 1:
                        y = value;
                        break;

                    case 2:
                        z = value;
                        break;

                    default:
                        throw new IndexOutOfRangeException();
                }
            }
        }

        public static Vec3 operator +(Vec3 a, Vec3 b)
        {
            return new Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        public static Vec3 operator -(Vec3 a, Vec3 b)
        {
            return new Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
        }

        public static Vec3 operator *(float f, Vec3 v)
        {
            return new Vec3(f * v.x, f * v.y, f * v.z);
        }

        public static Vec3 operator *(Vec3 v, float f)
        {
            return new Vec3(v.x * f, v.y * f, v.z * f);
        }

        public static Vec3 operator /(Vec3 v, float f)
        {
            return (1.0f / f) * v;
        }

        public static float Dot(Vec3 a, Vec3 b)
        {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }
    }

    public struct Vec4
    {
        public float x;
        public float y;
        public float z;
        public float w;

        public Vec4(float x, float y, float z, float w)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
        }

        public Vec4(float f)
            : this(f, f, f, f)
        { }

        public Vec4(Vec3 v, float w)
            : this(v.x, v.y, v.z, w)
        { }

        public static Vec4 operator +(Vec4 a, Vec4 b)
        {
            return new Vec4(a.x + b.x, a.y + b.y, a.x + b.z, a.w + b.w);
        }
    }

    public struct Mat3
    {
        public Vec3 col0;
        public Vec3 col1;
        public Vec3 col2;

        public Mat3(Vec3 col0, Vec3 col1, Vec3 col2)
        {
            this.col0 = col0;
            this.col1 = col1;
            this.col2 = col2;
        }

        public Mat3(float f)
            : this(
                new Vec3(f, 0.0f, 0.0f),
                new Vec3(0.0f, f, 0.0f),
                new Vec3(0.0f, 0.0f, f))
        { }

        public Mat3(
            float col0x, float col0y, float col0z,
            float col1x, float col1y, float col1z,
            float col2x, float col2y, float col2z)

            : this(
                new Vec3(col0x, col0y, col0z),
                new Vec3(col1x, col1y, col1z),
                new Vec3(col2x, col2y, col2z))
        { }


        public Vec3 this[int i]
        {
            get
            {
                switch (i)
                {
                    case 0:
                        return col0;

                    case 1:
                        return col1;

                    case 2:
                        return col2;

                    default:
                        throw new IndexOutOfRangeException();
                }
            }

            set
            {
                switch (i)
                {
                    case 0:
                        col0 = value;
                        break;

                    case 1:
                        col1 = value;
                        break;

                    case 2:
                        col2 = value;
                        break;

                    default:
                        throw new IndexOutOfRangeException();
                }
            }
        }

        public float this[int col, int row]
        {
            get
            {
                switch (col)
                {
                    case 0:
                        return col0[row];

                    case 1:
                        return col1[row];

                    case 2:
                        return col2[row];

                    default:
                        throw new IndexOutOfRangeException();
                }
            }

            set
            {
                switch (col)
                {
                    case 0:
                        col0[row] = value;
                        break;

                    case 1:
                        col1[row] = value;
                        break;

                    case 2:
                        col2[row] = value;
                        break;

                    default:
                        throw new IndexOutOfRangeException();
                }
            }
        }

        public static Vec3 operator *(Mat3 m, Vec3 v)
        {
            return new Vec3(
                m.col0.x * v.x + m.col1.x * v.y + m.col2.x * v.z,
                m.col0.y * v.x + m.col1.y * v.y + m.col2.y * v.z,
                m.col0.z * v.x + m.col1.z * v.y + m.col2.z * v.z);
        }

        public static Vec3 operator *(Vec3 v, Mat3 m)
        {
            return new Vec3(
                v.x * m.col0.x + v.y * m.col0.y + v.z * m.col0.z,
                v.x * m.col1.x + v.y * m.col1.y + v.z * m.col1.z,
                v.x * m.col2.x + v.y * m.col2.y + v.z * m.col2.z);
        }
    }

    public static class Svd
    {
        const int SVD_NUM_SWEEPS = 5;
        const float Tiny_Number = 1.0e-20f;

        public static float Abs(float x)
        {
            return Math.Abs(x);
        }

        public static float Sqrt(float x)
        {
            return (float)Math.Sqrt(x);
        }

        public static float Rsqrt(float x)
        {
            return 1.0f / (float)Math.Sqrt(x);
        }

        public static void GivensCoeffsSym(
            float app, float apq, float aqq, out float c, out float s)
        {
            if (apq == 0.0f)
            {
                c = 1.0f;
                s = 0.0f;
                return;
            }

            float tau = (aqq - app) / (2.0f * apq);
            float stt = Sqrt(1.0f + tau * tau);
            float tan = 1.0f / ((tau >= 0.0f) ? (tau + stt) : (tau - stt));
            c = Rsqrt(1.0f + tan * tan);
            s = tan * c;
        }

        public static void RotateXy(
            ref float x, ref float y, float c, float s)
        {
            float u = x; float v = y;
            x = c * u - s * v;
            y = s * u + c * v;
        }

        public static void RotateqXy(
            ref float x, ref float y, float a, float c, float s)
        {
            float cc = c * c;
            float ss = s * s;
            float mx = 2.0f * c * s * a;
            float u = x;
            float v = y;
            x = cc * u - mx + ss * v;
            y = ss * u + mx + cc * v;
        }

        public static void Rotate(ref Mat3 vtav, ref Mat3 v, int a, int b)
        {
            if (vtav[a, b] == 0.0f)
            {
                return;
            }

            float c, s;
            GivensCoeffsSym(vtav[a, a], vtav[a, b], vtav[b, b], out c, out s);

            {
                float x = vtav[a, a];
                float y = vtav[b, b];
                RotateqXy(ref x, ref y, vtav[a, b], c, s);
                vtav[a, a] = x;
                vtav[b, b] = y;
            }

            {
                float x = vtav[0, 3 - b];
                float y = vtav[1 - a, 2];
                RotateXy(ref x, ref y, c, s);
                vtav[0, 3 - b] = x;
                vtav[1 - a, 2] = y;
            }

            vtav[a, b] = 0.0f;

            {
                float x = v[0, a];
                float y = v[0, b];
                RotateXy(ref x, ref y, c, s);
                v[0, a] = x;
                v[0, b] = y;
            }

            {
                float x = v[1, a];
                float y = v[1, b];
                RotateXy(ref x, ref y, c, s);
                v[1, a] = x;
                v[1, b] = y;
            }

            {
                float x = v[2, a];
                float y = v[2, b];
                RotateXy(ref x, ref y, c, s);
                v[2, a] = x;
                v[2, b] = y;
            }
        }

        public static void SolveSym(Mat3 a, out Vec3 sigma, ref Mat3 v)
        {
            // assuming that A is symmetric: can optimize all operations for 
            // the upper right triagonal
            Mat3 vtav = a;

            // assuming V is identity: you can also pass a matrix the rotations
            // should be applied to
            // U is not computed
            for (int i = 0; i < SVD_NUM_SWEEPS; ++i)
            {
                Rotate(ref vtav, ref v, 0, 1);
                Rotate(ref vtav, ref v, 0, 2);
                Rotate(ref vtav, ref v, 1, 2);
            }

            sigma = new Vec3(vtav[0, 0], vtav[1, 1], vtav[2, 2]);
        }

        public static float Invdet(float x, float tol)
        {
            return (Abs(x) < tol || Abs(1.0f / x) < tol) ? 0.0f : (1.0f / x);
        }

        public static void Pseudoinverse(out Mat3 o, Vec3 sigma, Mat3 v)
        {
            float d0 = Invdet(sigma[0], Tiny_Number);
            float d1 = Invdet(sigma[1], Tiny_Number);
            float d2 = Invdet(sigma[2], Tiny_Number);

            o = new Mat3(
                v[0, 0] * d0 * v[0, 0] + v[0, 1] * d1 * v[0, 1] + v[0, 2] * d2 * v[0, 2],
                v[0, 0] * d0 * v[1, 0] + v[0, 1] * d1 * v[1, 1] + v[0, 2] * d2 * v[1, 2],
                v[0, 0] * d0 * v[2, 0] + v[0, 1] * d1 * v[2, 1] + v[0, 2] * d2 * v[2, 2],
                v[1, 0] * d0 * v[0, 0] + v[1, 1] * d1 * v[0, 1] + v[1, 2] * d2 * v[0, 2],
                v[1, 0] * d0 * v[1, 0] + v[1, 1] * d1 * v[1, 1] + v[1, 2] * d2 * v[1, 2],
                v[1, 0] * d0 * v[2, 0] + v[1, 1] * d1 * v[2, 1] + v[1, 2] * d2 * v[2, 2],
                v[2, 0] * d0 * v[0, 0] + v[2, 1] * d1 * v[0, 1] + v[2, 2] * d2 * v[0, 2],
                v[2, 0] * d0 * v[1, 0] + v[2, 1] * d1 * v[1, 1] + v[2, 2] * d2 * v[1, 2],
                v[2, 0] * d0 * v[2, 0] + v[2, 1] * d1 * v[2, 1] + v[2, 2] * d2 * v[2, 2]);
        }

        public static void SolveAtaAtb(Mat3 ata, Vec3 atb, out Vec3 x)
        {
            Mat3 V = new Mat3(1.0f);
            Vec3 sigma;

            SolveSym(ata, out sigma, ref V);

            // A = UEV^T; U = A / (E*V^T)
            Mat3 Vinv;
            Pseudoinverse(out Vinv, sigma, V);
            x = Vinv * atb;
        }

        public static Vec3 VmulSym(Mat3 a, Vec3 v)
        {
            return new Vec3(
                Vec3.Dot(a[0], v),
                a[0, 1] * v.x + a[1, 1] * v.y + a[1, 2] * v.z,
                a[0, 2] * v.x + a[1, 2] * v.y + a[2, 2] * v.z);
        }

        public static void MulAtaSym(out Mat3 o, Mat3 a)
        {
            o = new Mat3(0.0f);

            o[0, 0] = a[0, 0] * a[0, 0] + a[1, 0] * a[1, 0] + a[2, 0] * a[2, 0];
            o[0, 1] = a[0, 0] * a[0, 1] + a[1, 0] * a[1, 1] + a[2, 0] * a[2, 1];
            o[0, 2] = a[0, 0] * a[0, 2] + a[1, 0] * a[1, 2] + a[2, 0] * a[2, 2];
            o[1, 1] = a[0, 1] * a[0, 1] + a[1, 1] * a[1, 1] + a[2, 1] * a[2, 1];
            o[1, 2] = a[0, 1] * a[0, 2] + a[1, 1] * a[1, 2] + a[2, 1] * a[2, 2];
            o[2, 2] = a[0, 2] * a[0, 2] + a[1, 2] * a[1, 2] + a[2, 2] * a[2, 2];
        }

        public static void SolveAxB(
            Mat3 a, Vec3 b, out Mat3 ata, out Vec3 atb, out Vec3 x)
        {
            MulAtaSym(out ata, a);
            atb = b * a; // transpose(a) * b;
            SolveAtaAtb(ata, atb, out x);
        }
    }

    public static class Qef
    {
        public static void Add(
            Vec3 n, Vec3 p, ref Mat3 ata, ref Vec3 atb, ref Vec4 pointaccum)
        {
            ata[0, 0] += n.x * n.x;
            ata[0, 1] += n.x * n.y;
            ata[0, 2] += n.x * n.z;
            ata[1, 1] += n.y * n.y;
            ata[1, 2] += n.y * n.z;
            ata[2, 2] += n.z * n.z;

            float b = Vec3.Dot(p, n);
            atb += n * b;
            pointaccum += new Vec4(p, 1.0f);
        }

        public static float CalcError(Mat3 a, Vec3 x, Vec3 b)
        {
            Vec3 vtmp = b - Svd.VmulSym(a, x);
            return Vec3.Dot(vtmp, vtmp);
        }

        public static float Solve(
            Mat3 ata, Vec3 atb, Vec4 pointaccum,out Vec3 x)
        {
            Vec3 masspoint =
                new Vec3(pointaccum.x, pointaccum.y, pointaccum.z) /
                pointaccum.w;

            atb -= Svd.VmulSym(ata, masspoint);
            Svd.SolveAtaAtb(ata, atb, out x);
            float result = CalcError(ata, x, atb);

            x += masspoint;

            return result;
        }
    }
}
