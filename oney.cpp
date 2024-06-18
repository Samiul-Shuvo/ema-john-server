1. Bisection method
#include <iostream>
#include <iomanip>
#include <math.h>
#define f(x) x *x *x - x - 1
    using namespace std;
int main()
{
    float x0, x1, x, f0, f1, f, e;
    int step = 1;
    cout << setprecision(5) << fixed;
up:
    cout << "Enter first number: ";
    cin >> x0;
    cout << "Enter second number: ";
    cin >> x1;
    cout << "Enter tolerable error: ";
    cin >> e;
    f0 = f(x0);
    f1 = f(x1);
    if (f0 * f1 > 0.0)
    {
        cout << "Incorrect Initial Guesses." << endl;
        goto up;
    }
    cout << endl
         << "****************" << endl;
    cout << "Bisection Method BY C211066" << endl;
    cout << "****************" << endl;
    do
    {
        /* Bisecting Interval */
        x = (x0 + x1) / 2;
        f = f(x);
        cout << "Iteration-" << step << ":\t x = " << setw(10) << x << " and f(x) = " << setw(10) << f(x) << endl;
        if (f0 * f < 0)
        {
            x1 = x;
        }
        else
        {
            x0 = x;
        }
        step = step + 1;
    } while (fabs(f) > e);
    cout << endl
         << "Root is: " << x << endl;
    return 0;
}
2. False position method
#include <iostream>
#include <iomanip>
#include <math.h>
#define f(x) x *x *x - x - 1
    using namespace std;
int main()
{
    float x0, x1, x, f0, f1, f, e;
    int step = 1;
    cout << setprecision(5) << fixed;
up:
    cout << "Enter first number: ";
    cin >> x0;
    cout << "Enter second number: ";
    cin >> x1;
    cout << "Enter tolerable error: ";
    cin >> e;
    f0 = f(x0);
    f1 = f(x1);
    if (f0 * f1 > 0.0)
    {
        cout << "Incorrect Initial Guesses." << endl;
        goto up;
    }
    cout << endl
         << "*********************" << endl;
    cout << "False Position Method" << endl;
    cout << "*********************" << endl;
    do
    {
        x = x0 - (x0 - x1) * f0 / (f0 - f1);
        f = f(x);
        cout << "Iteration-" << step << ":\t x = " << setw(10) << x << " and f(x) = " << setw(10) << f(x) << endl;
        if (f0 * f < 0)
        {
            x1 = x;
            f1 = f;
        }
        else
        {
            x0 = x;
            f0 = f;
        }
        step = step + 1;
    } while (fabs(f) > e);
    cout << endl
         << "Root is: " << x << endl;
    return 0;
}
3. Newton - Raphson method
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#define f(x) x *x *x - x - 1
#define g(x) 3 * x *x - 1
    using namespace std;
int main()
{
    float x0, x1, f0, f1, g0, e;
    int step = 1, N;
    cout << setprecision(6) << fixed;
    cout << "Enter initial guess: ";
    cin >> x0;
    cout << "Enter tolerable error: ";
    cin >> e;
    cout << "Enter maximum iteration: ";
    cin >> N;
    cout << endl
         << "*********************" << endl;
    cout << "Newton Raphson Method" << endl;
    cout << "*********************" << endl;
    do
    {
        g0 = g(x0);
        f0 = f(x0);
        if (g0 == 0.0)
        {
            cout << "Mathematical Error.";
            exit(0);
        }
        x1 = x0 - f0 / g0;
        cout << "Iteration-" << step << ":\t x = " << setw(10) << x1 << " and f(x) = " << setw(10) << f(x1) << endl;
        x0 = x1;
        step = step + 1;
        if (step > N)
        {
            cout << "Not Convergent.";
            exit(0);
        }
        f1 = f(x1);
    } while (fabs(f1) > e);
    cout << endl
         << "Root is: " << x1;
    return 0;
}
4. Secand Method
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#define f(x) x *x - 4 * x - 10
    using namespace std;
int main()
{
    float x0, x1, x2, f0, f1, f2, e;
    int step = 1, N;
    cout << setprecision(6) << fixed;
    cout << "Enter first number: ";
    cin >> x0;
    cout << "Enter second number: ";
    cin >> x1;
    cout << "Enter tolerable error: ";
    cin >> e;
    cout << "Enter maximum iteration: ";
    cin >> N;
    cout << endl
         << "**************" << endl;
    cout << "Secant Method" << endl;
    cout << "**************" << endl;
    do
    {
        f0 = f(x0);
        f1 = f(x1);
        if (f0 == f1)
        {
            cout << "Mathematical Error.";
            exit(0);
        }
        x2 = x1 - (x1 - x0) * f1 / (f1 - f0);
        f2 = f(x2);
        cout << "Iteration-" << step << ":\t x2 = " << setw(10) << x2 << " and f(x2) = " << setw(10) << f(x2) << endl;
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
        step = step + 1;
        if (step > N)
        {
            cout << "Not Convergent.";
            exit(0);
        }
    } while (fabs(f2) > e);
    cout << endl
         << "Root is: " << x2;
    return 0;
}
5. Fixed point method
#include <iostream>
#include <cmath>
#include <iomanip>
    using namespace std;
double g(double x)
{
    return cbrt(x + 1);
}
int main()
{
    double x0, x1, tolerance;
    int max_iterations, iteration = 0;
    cout << "Enter initial number: ";
    cin >> x0;
    cout << "Enter error: ";
    cin >> tolerance;
    cout << "Enter maximum number of iterations: ";
    cin >> max_iterations;
    cout << fixed << setprecision(6);
    do
    {
        x1 = g(x0);
        iteration++;
        cout << "Iteration " << iteration << ": x = " << x1 << endl;
        if (fabs(x1 - x0) < tolerance)
        {
            cout << "Converged to root: " << x1 << endl;
            return 0;
        }
        x0 = x1;
    } while (iteration < max_iterations);
    cout << "Did not converge within the maximum number of iterations." << endl;
    return 0;
}
1. 6. Newton Forward Interpolation
#include <iostream>
#include <cmath>
#include <iomanip>
    double
    f(double x)
{
    return x * x * x - x - 1;
}
double f_prime(double x)
{
    return 3 * x * x - 1;
}
int main()
{
    double x0, x1, tolerance;
    int max_iterations, iteration = 0;
    std::cout << "Enter initial number: ";
    std::cin >> x0;
    std::cout << "Enter error: ";
    std::cin >> tolerance;
    std::cout << "Enter maximum number of iterations: ";
    std::cin >> max_iterations;
    std::cout << std::fixed << std::setprecision(6);
    do
    {
        x1 = x0 - f(x0) / f_prime(x0);
        iteration++;
        std::cout << "Iteration " << iteration << ": x = " << x1 << std::endl;
        if (std::fabs(x1 - x0) < tolerance)
        {
            std::cout << "Converged to root: " << x1 << std::endl;
            return 0;
        }
        x0 = x1;
    } while (iteration < max_iterations);
    std::cout << "Did not converge within the maximum number of iterations." << std::endl;
    return 0;
}
7. Newton Backward Interpolation
#include <iostream>
    using namespace std;
int main()
{
    float x[20], y[20][20];
    int i, j, n;
    cout << "Enter number of data? " << endl;
    cin >> n;
    cout << "Enter data: " << endl;
    for (i = 0; i < n; i++)
    {
        cout << "x[" << i << "] = ";
        cin >> x[i];
        cout << "y[" << i << "] = ";
        cin >> y[i][0];
    }
    for (i = 1; i < n; i++)
    {
        for (j = n - 1; j > i - 1; j--)
        {
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
        }
    }
    cout << endl
         << "BACKWARD DIFFERENCE TABLE" << endl;
    for (i = 0; i < n; i++)
    {
        cout << x[i];
        for (j = 0; j <= i; j++)
        {
            cout << "\t" << y[i][j];
        }
        cout << endl;
    }
    return 0;
}
8. Lagranges interpolation
#include <iostream>
#include <conio.h>
    using namespace std;
int main()
{
    float x[100], y[100], xp, yp = 0, p;
    int i, j, n;
    cout << "Enter number of data: ";
    cin >> n;
    cout << "Enter data:" << endl;
    for (i = 1; i <= n; i++)
    {
        cout << "x[" << i << "] = ";
        cin >> x[i];
        cout << "y[" << i << "] = ";
        cin >> y[i];
    }
    cout << "Enter interpolation point: ";
    cin >> xp;
    for (i = 1; i <= n; i++)
    {
        p = 1;
        for (j = 1; j <= n; j++)
        {
            if (i != j)
            {
                p = p * (xp - x[j]) / (x[i] - x[j]);
            }
        }
        yp = yp + p * y[i];
    }
    cout << endl
         << "Interpolated value at " << xp << " is " << yp;
    return 0;
}
9. Newtons Divided Differnece
#include <iostream>
    using namespace std;
int main()
{
    float x[20], y[20][20];
    int i, j, n;
    cout << "Enter number of data? " << endl;
    cin >> n;
    cout << "Enter data: " << endl;
    for (i = 0; i < n; i++)
    {
        cout << "x[" << i << "] = ";
        cin >> x[i];
        cout << "y[" << i << "] = ";
        cin >> y[i][0];
    }
    for (i = 1; i < n; i++)
    {
        for (j = 0; j < n - i; j++)
        {
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }
    }
    cout << endl
         << "FORWARD DIFFERENCE TABLE" << endl;
    for (i = 0; i < n; i++)
    {
        cout << x[i];
        for (j = 0; j < n - i; j++)
        {
            cout << "\t" << y[i][j];
        }
        cout << endl;
    }
    return 0;
}
10. Matrix Inversion method
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#define SIZE 10
    using namespace std;
int main()
{
    float a[SIZE][SIZE], x[SIZE], ratio;
    int i, j, k, n;
    cout << setprecision(3) << fixed;
    cout << "Enter order of matrix: ";
    cin >> n;
    cout << "Enter coefficients of Matrix: " << endl;
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            cout << "a[" << i << "]" << j << "]= ";
            cin >> a[i][j];
        }
    }
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            if (i == j)
            {
                a[i][j + n] = 1;
            }
            else
            {
                a[i][j + n] = 0;
            }
        }
    }
    for (i = 1; i <= n; i++)
    {
        if (a[i][i] == 0.0)
        {
            cout << "Mathematical Error!";
            exit(0);
        }
        for (j = 1; j <= n; j++)
        {
            if (i != j)
            {
                ratio = a[j][i] / a[i][i];
                for (k = 1; k <= 2 * n; k++)
                {
                    a[j][k] = a[j][k] - ratio * a[i][k];
                }
            }
        }
    }
    for (i = 1; i <= n; i++)
    {
        for (j = n + 1; j <= 2 * n; j++)
        {
            a[i][j] = a[i][j] / a[i][i];
        }
    }
    cout << endl
         << "Inverse Matrix is:" << endl;
    for (i = 1; i <= n; i++)
    {
        for (j = n + 1; j <= 2 * n; j++)
        {
            cout << a[i][j] << "\t";
        }
        cout << endl;
    }
    return (0);
}
11. Gauss Elimination method
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#define SIZE 10
    using namespace std;
int main()
{
    float a[SIZE][SIZE], x[SIZE], ratio;
    int i, j, k, n;
    cout << setprecision(3) << fixed;
    cout << "Enter number of unknowns: ";
    cin >> n;
    cout << "Enter Coefficients of Augmented Matrix: " << endl;
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n + 1; j++)
        {
            cout << "a[" << i << "]" << j << "]= ";
            cin >> a[i][j];
        }
    }
    for (i = 1; i <= n - 1; i++)
    {
        if (a[i][i] == 0.0)
        {
            cout << "Mathematical Error!";
            exit(0);
        }
        for (j = i + 1; j <= n; j++)
        {
            ratio = a[j][i] / a[i][i];
            for (k = 1; k <= n + 1; k++)
            {
                a[j][k] = a[j][k] - ratio * a[i][k];
            }
        }
    }
    x[n] = a[n][n + 1] / a[n][n];
    for (i = n - 1; i >= 1; i--)
    {
        x[i] = a[i][n + 1];
        for (j = i + 1; j <= n; j++)
        {
            x[i] = x[i] - a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }
    cout << endl
         << "Solution: " << endl;
    for (i = 1; i <= n; i++)
    {
        cout << "x[" << i << "] = " << x[i] << endl;
    }
    return (0);
}
12. Gauss Jordan method
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#define SIZE 10
    using namespace std;
int main()
{
    float a[SIZE][SIZE], x[SIZE], ratio;
    int i, j, k, n;
    cout << setprecision(3) << fixed;
    cout << "Enter number of unknowns: ";
    cin >> n;
    cout << "Enter Coefficients of Augmented Matrix: " << endl;
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n + 1; j++)
        {
            cout << "a[" << i << "]" << j << "]= ";
            cin >> a[i][j];
        }
    }
    for (i = 1; i <= n; i++)
    {
        if (a[i][i] == 0.0)
        {
            cout << "Mathematical Error!";
            exit(0);
        }
        for (j = 1; j <= n; j++)
        {
            if (i != j)
            {
                ratio = a[j][i] / a[i][i];
                for (k = 1; k <= n + 1; k++)
                {
                    a[j][k] = a[j][k] - ratio * a[i][k];
                }
            }
        }
    }
    for (i = 1; i <= n; i++)
    {
        x[i] = a[i][n + 1] / a[i][i];
    }
    cout << endl
         << "Solution: " << endl;
    for (i = 1; i <= n; i++)
    {
        cout << "x[" << i << "] = " << x[i] << endl;
    }
    return (0);
}
13. Dolittle LU method
#include <iostream>
    using namespace std;
void LUdecomposition(float a[10][10], float l[10][10], float u[10][10], int n)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j < i)
                l[j][i] = 0;
            else
            {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++)
                {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++)
        {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else
            {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++)
                {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
}
int main()
{
    float a[10][10], l[10][10], u[10][10];
    int n = 0, i = 0, j = 0;
    cout << "Enter size of square matrix : " << endl;
    cin >> n;
    cout << "Enter matrix values: " << endl;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            cin >> a[i][j];
    LUdecomposition(a, l, u, n);
    cout << "L Decomposition is as follows..." << endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            cout << l[i][j] << " ";
        }
        cout << endl;
    }
    cout << "U Decomposition is as follows..." << endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            cout << u[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}
14. Jacobi's method
#include <iostream>
#include <iomanip>
#include <math.h>
#define f1(x, y, z) (17 - y + 2 * z) / 20
#define f2(x, y, z) (-18 - 3 * x + z) / 20
#define f3(x, y, z) (25 - 2 * x + 3 * y) / 20
    using namespace std;
int main()
{
    float x0 = 0, y0 = 0, z0 = 0, x1, y1, z1, e1, e2, e3, e;
    int step = 1;
    cout << setprecision(6) << fixed;
    cout << "Enter tolerable error: ";
    cin >> e;
    cout << endl
         << "Count\tx\t\ty\t\tz" << endl;
    do
    {
        x1 = f1(x0, y0, z0);
        y1 = f2(x0, y0, z0);
        z1 = f3(x0, y0, z0);
        cout << step << "\t" << x1 << "\t" << y1 << "\t" << z1 << endl;
        e1 = fabs(x0 - x1);
        e2 = fabs(y0 - y1);
        e3 = fabs(z0 - z1);
        step++;
        x0 = x1;
        y0 = y1;
        z0 = z1;
    } while (e1 > e && e2 > e && e3 > e);
    cout << endl
         << "Solution: x = " << x1 << ", y = " << y1 << " and z = " << z1;
    return 0;
}

15. Gauss - Seidel method

#include <iostream>
#include <iomanip>
#include <math.h>
#define f1(x, y, z) (17 - y + 2 * z) / 20
#define f2(x, y, z) (-18 - 3 * x + z) / 20
#define f3(x, y, z) (25 - 2 * x + 3 * y) / 20
    using namespace std;
int main()
{
    float x0 = 0, y0 = 0, z0 = 0, x1, y1, z1, e1, e2, e3, e;
    int step = 1;
    cout << setprecision(6) << fixed;
    cout << "Enter tolerable error: ";
    cin >> e;
    cout << endl
         << "Count\tx\t\ty\t\tz" << endl;
    do
    {
        x1 = f1(x0, y0, z0);
        y1 = f2(x1, y0, z0);
        z1 = f3(x1, y1, z0);
        cout << step << "\t" << x1 << "\t" << y1 << "\t" << z1 << endl;
        e1 = fabs(x0 - x1);
        e2 = fabs(y0 - y1);
        e3 = fabs(z0 - z1);
        step++;
        x0 = x1;
        y0 = y1;
        z0 = z1;
    } while (e1 > e && e2 > e && e3 > e);
    cout << endl
         << "Solution: x = " << x1 << ", y = " << y1 << " and z = " << z1;
    return 0;
}
16. Least Square method
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
    using namespace std;
void leastSquaresQuadratic(const vector<double> &x, const vector<double> &y, double &a, double &b, double &c)
{
    int n = x.size();
    double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_x3 = 0, sum_x4 = 0;
    double sum_xy = 0, sum_x2y = 0;
    for (int i = 0; i < n; ++i)
    {
        double xi = x[i];
        double yi = y[i];
        double xi2 = xi * xi;
        double xi3 = xi2 * xi;
        double xi4 = xi3 * xi;
        sum_x += xi;
        sum_y += yi;
        sum_x2 += xi2;
        sum_x3 += xi3;
        sum_x4 += xi4;
        sum_xy += xi * yi;
        sum_x2y += xi2 * yi;
    }
    double D = n * (sum_x2 * sum_x4 - sum_x3 * sum_x3) - sum_x * (sum_x * sum_x4 - sum_x2 * sum_x3) + sum_x2 * (sum_x * sum_x3 - sum_x2 * sum_x2);
    double Da = sum_y * (sum_x2 * sum_x4 - sum_x3 * sum_x3) - sum_x * (sum_xy * sum_x4 - sum_x2y * sum_x3) +
                sum_x2 * (sum_xy * sum_x3 - sum_x2y * sum_x2);
    double Db = n * (sum_xy * sum_x4 - sum_x2y * sum_x3) - sum_y * (sum_x * sum_x4 - sum_x2 * sum_x3) +
                sum_x2 * (sum_x * sum_x2y - sum_x2 * sum_xy);
    double Dc = n * (sum_x2 * sum_x2y - sum_x3 * sum_xy) - sum_x * (sum_x * sum_x2y - sum_x2 * sum_xy) +
                sum_y * (sum_x * sum_x3 - sum_x2 * sum_x2);
    a = Da / D;
    b = Db / D;
    c = Dc / D;
}
int main()
{
    int n;
    cout << "Enter the number of data points: ";
    cin >> n;
    vector<double> x(n), y(n);
    cout << "Enter the x and y values:" << endl;
    for (int i = 0; i < n; ++i)
    {
        cout << "x[" << i << "]: ";
        cin >> x[i];
        cout << "y[" << i << "]: ";
        cin >> y[i];
    }
    double a, b, c;
    leastSquaresQuadratic(x, y, a, b, c);
    cout << fixed << setprecision(6);
    cout << "The best fitting quadratic curve is y = " << a << "x^2 + " << b << "x + " << c << endl;
    cout << "Fitted values:" << endl;
    for (int i = 0; i < n; ++i)
    {
        double y_fit = a * x[i] * x[i] + b * x[i] + c;
        cout << "x = " << x[i] << ", y_fit = " << y_fit << endl;
    }
    return 0;
}
17. Derivatives using Newton's Forward interpolation
#include <iostream>
    using namespace std;
int main()
{
    float x[20], y[20][20];
    int i, j, n;
    cout << "Enter number of data? " << endl;
    cin >> n;
    cout << "Enter data: " << endl;
    for (i = 0; i < n; i++)
    {
        cout << "x[" << i << "] = ";
        cin >> x[i];
        cout << "y[" << i << "] = ";
        cin >> y[i][0];
    }
    for (i = 1; i < n; i++)
    {
        for (j = 0; j < n - i; j++)
        {
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }
    }
    cout << endl
         << "FORWARD DIFFERENCE TABLE" << endl;
    for (i = 0; i < n; i++)
    {
        cout << x[i];
        for (j = 0; j < n - i; j++)
        {
            cout << "\t" << y[i][j];
        }
        cout << endl;
    }
    return 0;
}
18. Derivaties using Backward difference interpolation.
#include <iostream>
#include <vector>
#include <iomanip>
    using namespace std;
void buildBackwardDifferenceTable(const vector<double> &x, const vector<double> &y, vector<vector<double>> &diffTable)
{
    int n = x.size();
    for (int i = 0; i < n; ++i)
    {
        diffTable[i][0] = y[i];
    }
    for (int j = 1; j < n; ++j)
    {
        for (int i = n - 1; i >= j; --i)
        {
            diffTable[i][j] = diffTable[i][j - 1] - diffTable[i - 1][j - 1];
        }
    }
}
double backwardDifferenceDerivative(const vector<double> &x, const vector<vector<double>> &diffTable, double h,
                                    double value)
{
    int n = x.size();
    int index = n - 1;
    for (int i = 0; i < n; ++i)
    {
        if (x[i] > value)
        {
            index = i - 1;
            break;
        }
    }
    double u = (value - x[index]) / h;
    double derivative = diffTable[index][1] / h; // First term (first backward difference)
    for (int i = 2; i < n; ++i)
    {
        double term = diffTable[index][i];
        for (int j = 0; j < i - 1; ++j)
        {
            term *= (u + j);
        }
        derivative += term / (h * i);
    }
    return derivative;
}
int main()
{
    int n;
    cout << "Enter the number of data points: ";
    cin >> n;
    vector<double> x(n), y(n);
    cout << "Enter the x and y values:" << endl;
    for (int i = 0; i < n; ++i)
    {
        cout << "x[" << i << "]: ";
        cin >> x[i];
        cout << "y[" << i << "]: ";
        cin >> y[i];
    }
    vector<vector<double>> diffTable(n, vector<double>(n, 0.0));
    buildBackwardDifferenceTable(x, y, diffTable);
    double value, h;
    cout << "Enter the point at which the derivative is to be computed: ";
    cin >> value;
    cout << "Enter the step size h: ";
    cin >> h;
    double derivative = backwardDifferenceDerivative(x, diffTable, h, value);
    cout << fixed << setprecision(6);
    cout << "The derivative at x = " << value << " using backward difference interpolation is " << derivative << endl;
    return 0;
}
19. Derivatives using Central difference interpolation
#include <iostream>
#include <cmath>
    using namespace std;
double function(double x)
{
    return x * x; // Example function: f(x) = x^2
}
double centralDifference(double (*f)(double), double x, double h)
{
    return (f(x + h) - f(x - h)) / (2 * h);
}
int main()
{
    double x, h;
    cout << "Enter the point x at which the derivative is to be computed: ";
    cin >> x;
    cout << "Enter the step size h: ";
    cin >> h;
    double derivative = centralDifference(function, x, h);
    cout << "The central difference derivative at x = " << x << " is " << derivative << endl;
    return 0;
}
20. Taylor's series method
#include <iostream>
#include <cmath>
    using namespace std;
double function(double x)
{
    return exp(x);
}
double nthDerivative(double (*f)(double), double x, int n)
{
    if (n == 0)
    {
        return f(x);
    }
    else
    {
        return nthDerivative(f, x + 0.0001, n - 1) - nthDerivative(f, x - 0.0001, n - 1);
    }
}
int factorial(int n);
double taylorSeries(double (*f)(double), double x, int n)
{
    double result = 0;
    for (int i = 0; i < n; ++i)
    {
        result += nthDerivative(f, x, i) / factorial(i) * pow(x, i);
    }
    return result;
}
int factorial(int n)
{
    if (n <= 1)
    {
        return 1;
    }
    return n * factorial(n - 1);
}
int main()
{
    double x;
    int n;
    cout << "Enter the point x at which the Taylor series expansion is to be computed: ";
    cin >> x;
    cout << "Enter the number of terms in the series: ";
    cin >> n;
    double result = taylorSeries(function, x, n);
    cout << "The value of the Taylor series expansion at x = " << x << " using " << n << " terms is " << result << endl;
    return 0;
}
21. Euler's method
#include <iostream>
#define f(x, y) x + y
    using namespace std;
int main()
{
    float x0, y0, xn, h, yn, slope;
    int i, n;
    cout << "Enter Initial Condition" << endl;
    cout << "x0 = ";
    cin >> x0;
    cout << "y0 = ";
    cin >> y0;
    cout << "Enter calculation point xn = ";
    cin >> xn;
    cout << "Enter number of steps: ";
    cin >> n;
    h = (xn - x0) / n;
    cout << "\nx0\ty0\tslope\tyn\n";
    cout << "------------------------------\n";
    for (i = 0; i < n; i++)
    {
        slope = f(x0, y0);
        yn = y0 + h * slope;
        cout << x0 << "\t" << y0 << "\t" << slope << "\t" << yn << endl;
        y0 = yn;
        x0 = x0 + h;
    }
    cout << "\nValue of y at x = " << xn << " is " << yn;
    return 0;
}
22. Runge - Kutta method
#include <iostream>
#define f(x, y) (y * y - x * x) / (y * y + x * x)
    using namespace std;
#define f(x, y) (y * y - x * x) / (y * y + x * x)
using namespace std;
int main()
{
    float x0, y0, xn, h, yn, k1, k2, k3, k4, k;
    int i, n;
    cout << "Enter Initial Condition" << endl;
    cout << "x0 = ";
    cin >> x0;
    cout << "y0 = ";
    cin >> y0;
    cout << "Enter calculation point xn = ";
    cin >> xn;
    cout << "Enter number of steps: ";
    cin >> n;
    h = (xn - x0) / n;
    cout << "\nx0\ty0\tyn\n";
    cout << "------------------\n";
    for (i = 0; i < n; i++)
    {
        k1 = h * (f(x0, y0));
        k2 = h * (f((x0 + h / 2), (y0 + k1 / 2)));
        k3 = h * (f((x0 + h / 2), (y0 + k2 / 2)));
        k4 = h * (f((x0 + h), (y0 + k3)));
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        yn = y0 + k;
        cout << x0 << "\t" << y0 << "\t" << yn << endl;
        x0 = x0 + h;
        y0 = yn;
    }
    cout << "\nValue of y at x = " << xn << " is " << yn;
    return 0;
}