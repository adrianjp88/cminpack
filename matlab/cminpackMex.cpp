#include <iostream>
#include "mex.h"
#include <stdlib.h>
#include <vector>

#define real __cminpack_real__
#define __cminpack_float__
#include "..\cminpack.h"

#define GAUSS_1D 0
#define GAUSS_2D 1
#define GAUSS_2D_ELLIPTIC 2
#define GAUSS_2D_ROTATED 3
#define CAUCHY_2D_ELLIPTIC 4
#define LINEAR_1D 5
#define FLETCHER_POWELL 6
#define BROWN_DENNIS 7
#define	RAMSEY_FIXED_P 8
#define	RAMSEY_VAR_P 9


typedef struct  {
    int size;
    real * data;
    real * x_data;
} fcndata_t;

int gauss1d(void *p,
    int size,
    int n_parameters,
    real const * parameters,
    real * function,
    real * jacobian,
    int ldfjac,
    int iflag)
{
    const real *data = ((fcndata_t*)p)->data;

    for (int x = 0; x < size; x++)
    {
        real const arg
            = ((x - parameters[1])*(x - parameters[1]))
            / (2.f * parameters[2] * parameters[2]);

        real ex = exp(-arg);

        if (iflag != 2)
        {
            function[x] = data[x] - (parameters[0] * ex + parameters[3]);
        }
        else
        {
            jacobian[x + size * 0] = -ex;
            jacobian[x + size * 1] = -(parameters[0] * (x - parameters[1])*ex) / (parameters[2] * parameters[2]);
            jacobian[x + size * 2] = -(parameters[0] * (x - parameters[1])*(x - parameters[1])*ex) / (parameters[2] * parameters[2] * parameters[2]);
            jacobian[x + size * 3] = -1.f;
        }
    }
    return 0;
}

int gauss2d_elliptic(void *p,
    int size,
    int n_parameters,
    real const * parameters,
    real * function,
    real * jacobian,
    int ldfjac,
    int iflag)
{
    const real *data = ((fcndata_t*)p)->data;
    int const size_x = (int)std::sqrt(size);
    int const size_y = (int)std::sqrt(size);

    if (iflag != 2)
    {
        for (int iy = 0; iy < size_y; iy++)
        {
            for (int ix = 0; ix < size_x; ix++)
            {
                int const pixel_index = iy * size_x + ix;

                real const argx
                    = exp(-0.5f * ((ix - parameters[1]) / parameters[3]) * ((ix - parameters[1]) / parameters[3]));
                real const argy
                    = exp(-0.5f * ((iy - parameters[2]) / parameters[4]) * ((iy - parameters[2]) / parameters[4]));

                function[pixel_index] = data[pixel_index] - (parameters[0] * argx * argy + parameters[5]);
            }
        }
    }
    else
    {
        for (int iy = 0; iy < size_y; iy++)
        {
            for (int ix = 0; ix < size_x; ix++)
            {
                real argx = ((ix - parameters[1])*(ix - parameters[1])) / (2 * parameters[3] * parameters[3]);
                real argy = ((iy - parameters[2])*(iy - parameters[2])) / (2 * parameters[4] * parameters[4]);
                real ex = exp(-argx) * exp(-argy);

                int const pixel_index = iy*size_x + ix;

                jacobian[pixel_index + size * 0] = -ex;
                jacobian[pixel_index + size * 1] = -(parameters[0] * (ix - parameters[1])*ex) / (parameters[3] * parameters[3]);
                jacobian[pixel_index + size * 2] = -(parameters[0] * (iy - parameters[2])*ex) / (parameters[4] * parameters[4]);
                jacobian[pixel_index + size * 3] = -(parameters[0] * (ix - parameters[1])*(ix - parameters[1])*ex) / (parameters[3] * parameters[3] * parameters[3]);
                jacobian[pixel_index + size * 4] = -(parameters[0] * (iy - parameters[2])*(iy - parameters[2])*ex) / (parameters[4] * parameters[4] * parameters[4]);
                jacobian[pixel_index + size * 5] = -1;
            }
        }
    }
    return 0;
}

int gauss2d(void *p,
    int size,
    int n_parameters,
    real const * parameters,
    real * function,
    real * jacobian,
    int ldfjac,
    int iflag)
{
    const real *data = ((fcndata_t*)p)->data;
    int const size_x = (int)std::sqrt(size);
    int const size_y = (int)std::sqrt(size);

    if (iflag != 2)
    {
        for (int iy = 0; iy < size_y; iy++)
        {
            for (int ix = 0; ix < size_x; ix++)
            {
                int const pixel_index = iy * size_x + ix;

                real const argx
                    = exp(-0.5f * ((ix - parameters[1]) / parameters[3]) * ((ix - parameters[1]) / parameters[3]));
                real const argy
                    = exp(-0.5f * ((iy - parameters[2]) / parameters[3]) * ((iy - parameters[2]) / parameters[3]));

                function[pixel_index] = data[pixel_index] - (parameters[0] * argx * argy + parameters[4]);
            }
        }
    }
    else
    {
        for (int iy = 0; iy < size_y; iy++)
        {
            for (int ix = 0; ix < size_x; ix++)
            {
                real argx = ((ix - parameters[1])*(ix - parameters[1])) / (2 * parameters[3] * parameters[3]);
                real argy = ((iy - parameters[2])*(iy - parameters[2])) / (2 * parameters[3] * parameters[3]);
                real ex = exp(-argx) * exp(-argy);

                int const pixel_index = iy*size_x + ix;

                jacobian[pixel_index + size * 0] = -ex;
                jacobian[pixel_index + size * 1] = -(parameters[0] * (ix - parameters[1])*ex) / (parameters[3] * parameters[3]);
                jacobian[pixel_index + size * 2] = -(parameters[0] * (iy - parameters[2])*ex) / (parameters[3] * parameters[3]);
                jacobian[pixel_index + size * 3]
                    = -(parameters[0] * ex
                    * ((ix - parameters[1])*(ix - parameters[1])
                    + (iy - parameters[2])*(iy - parameters[2])))
                    / (parameters[3] * parameters[3] * parameters[3]);
                jacobian[pixel_index + size * 4] = -1;
            }
        }
    }
    return 0;
}

int calculate_ramsey_fixed_p(
    void * model,
    int size,
    int n_parameters,
    real const * parameters,
    real * function,
    real * jacobian,
    int ldfjac,
    int iflag)
{
    const real * data = ((fcndata_t*)model)->data;
    const real * x_data = ((fcndata_t*)model)->x_data;
    real const * p = parameters;

    real const pi = 3.14159f;

    //////////////////////////// values //////////////////////////////

    // parameters: [A1 A2 c f1 f2 t2star x1 x2] exp(-(point_index./t2star)^1)*(A1*cos(2*pi*f1*(point_index - x1)) + A2*cos(2*pi*f2*(point_index-x2))) + c

    for (int i = 0; i < size; i++)
    {
        real x;

        if (!x_data)
            x = real(i);
        else
            x = x_data[i];

        real const t2arg = x / p[5];
        real const ex = exp(-t2arg);
        real const phasearg1 = 2.0f * pi*p[3] * (x - p[6]);
        real const phasearg2 = 2.0f * pi*p[4] * (x - p[7]);
        real const cos1 = cos(phasearg1);
        real const sin1 = sin(phasearg1);
        real const cos2 = cos(phasearg2);
        real const sin2 = sin(phasearg2);

        if (iflag != 2)
        {
            function[i] = ex*(p[0] * cos1 + p[1] * cos2) + p[2] - data[i]; // formula calculating fit model values
        }
        else
        {
            /////////////////////////// derivatives ///////////////////////////
            real * current_jacobian = jacobian + i;

            current_jacobian[0 * size] = ex*cos1; // formula calculating derivative values with respect to parameters[0]
            current_jacobian[1 * size] = ex*cos2;
            current_jacobian[2 * size] = 1.0f;
            current_jacobian[3 * size] = -p[0] * 2.0f * pi*(x - p[6])*ex*sin1;
            current_jacobian[4 * size] = -p[1] * 2.0f * pi*(x - p[7])*ex*sin2;
            current_jacobian[5 * size] = 1.0f / (p[5] * p[5])*x*ex*(p[0] * cos1 + p[1] * cos2);
            current_jacobian[6 * size] = p[0] * 2.0f * pi*p[3] * sin1*ex;
            current_jacobian[7 * size] = p[1] * 2.0f * pi*p[4] * sin2*ex;
        }
    }
    return 0;
}

int calculate_ramsey_var_p(
    void * model,
    int size,
    int n_parameters,
    real const * parameters,
    real * function,
    real * jacobian,
    int ldfjac,
    int iflag)
{
    const real * data = ((fcndata_t*)model)->data;
    const real * x_data = ((fcndata_t*)model)->x_data;
    real const * p = parameters;

    real const pi = 3.14159f;

    ///////////////////////////// values //////////////////////////////

    // parameters: [A1 A2 c f1 f2 p t2star x1 x2] exp(-(x./t2star)^p)*(A1*cos(2*pi*f1*(x - x1)) + A2*cos(2*pi*f2*(x-x2))) + c
    for (int i = 0; i < size; i++)
    {
        real x;

        if (!x_data)
            x = real(i);
        else
            x = x_data[i];

        real const pi = 3.14159f;
        real const t2arg = pow(x / p[6], p[5]);
        real const ex = exp(-t2arg);
        real const phasearg1 = 2.f * pi*p[3] * (x - p[7]);
        real const phasearg2 = 2.f * pi*p[4] * (x - p[8]);
        real const cos1 = cos(phasearg1);
        real const sin1 = sin(phasearg1);
        real const cos2 = cos(phasearg2);
        real const sin2 = sin(phasearg2);
        //real const xmin = x/p[6] - 1;
        //real const log = xmin - xmin*xmin/2.f + xmin*xmin*xmin/3.f - xmin*xmin*xmin*xmin/4.f;

        if (iflag != 2)
        {
            function[i] = ex*(p[0] * cos1 + p[1] * cos2) + p[2] - data[i]; // formula calculating fit model values
        }
        else
        {
            /////////////////////////// derivatives ///////////////////////////
            real * current_jacobian = jacobian + i;
            current_jacobian[0 * size] = ex*cos1;
            current_jacobian[1 * size] = ex*cos2;
            current_jacobian[2 * size] = 1.f;
            current_jacobian[3 * size] = -p[0] * 2.f * pi*(x - p[7])*ex*sin1;
            current_jacobian[4 * size] = -p[1] * 2.f * pi*(x - p[8])*ex*sin2;
            current_jacobian[5 * size] = -log(x / p[6] + 0.000001f)*ex*t2arg*(p[0] * cos1 + p[1] * cos2);
            current_jacobian[6 * size] = p[5] * 1.f / (p[6] * p[6])*x*ex*pow(x / p[6], p[5] - 1.f)*(p[0] * cos1 + p[1] * cos2);
            current_jacobian[7 * size] = p[0] * 2.f * pi*p[3] * sin1*ex;
            current_jacobian[8 * size] = p[1] * 2.f * pi*p[4] * sin2*ex;
        }
    }
    return 0;
}

int calculate_fletcher_powell(
    void * model,
    int size,
    int n_parameters,
    real const * parameters,
    real * function,
    real * jacobian,
    int ldfjac,
    int iflag)
{
    real const * p = parameters;

    real const pi = 3.14159f;

    real theta = 0.f;

    if (p[0] > .0f)
        theta = .5f * atan(p[1] / p[0]) / pi;
    else if (p[0] < 0.)
        theta = .5f * atan(p[1] / p[0]) / pi + .5f;

    real const arg = p[0] * p[0] + p[1] * p[1];

    if (iflag != 2)
    {
        function[0] = 10.0f * (p[2] - 10.0f * theta);
        function[1] = 10.0f * (std::sqrt(arg) - 1.0f);
        function[2] = p[2];
    }
    else
    {
        // derivatives with respect to p[0]
        jacobian[0 * size + 0] = 100.0f * 1.0f / (2.0f*pi) * p[1] / arg;
        jacobian[0 * size + 1] = 10.0f * p[0] / std::sqrt(arg);
        jacobian[0 * size + 2] = 0.0f;

        // jacobian with respect to p[1]
        jacobian[1 * size + 0] = -100.0f * 1.0f / (2.0f*pi) * p[0] / (arg);
        jacobian[1 * size + 1] = 10.0f * p[1] / std::sqrt(arg);
        jacobian[1 * size + 2] = 0.0f;

        // jacobian with respect to p[2]
        jacobian[2 * size + 0] = 10.0f;
        jacobian[2 * size + 1] = 0.0f;
        jacobian[2 * size + 2] = 1.0f;
    }

    return 0;
}

int calculate_brown_dennis(
    void * model,
    int size,
    int n_parameters,
    real const * parameters,
    real * function,
    real * jacobian,
    int ldfjac,
    int iflag)
{
    real const * p = parameters;

    for (int point_index = 0; point_index < size; point_index++)
    {
        real const t = static_cast<real>(point_index) / 5.0f;

        real const arg1 = p[0] + p[1] * t - std::exp(t);
        real const arg2 = p[2] + p[3] * std::sin(t) - std::cos(t);

        if (iflag != 2)
        {
            function[point_index] = arg1*arg1 + arg2*arg2;
        }
        else
        {
            real * current_jacobian = jacobian + point_index;

            current_jacobian[0 * size] = 2.0f * arg1;
            current_jacobian[1 * size] = 2.0f * t * arg1;
            current_jacobian[2 * size] = 2.0f * arg2;
            current_jacobian[3 * size] = 2.0f * std::sin(t) * arg2;
        }
    }
    return 0;
}

void mexFunction(
    int          nlhs,
    mxArray      *plhs[],
    int          nrhs,
    mxArray const *prhs[])
{
    int expected_nrhs = 9;
    int expected_nlhs = 3;

    bool wrong_nrhs = false;
    bool wrong_nlhs = false;

    if (nrhs != expected_nrhs)
    {
        wrong_nrhs = true;
    }
    else if (nlhs != expected_nlhs)
    {
        wrong_nlhs = true;
    }

    if (wrong_nrhs || wrong_nlhs)
    {
        if (nrhs != expected_nrhs)
        {
            char s1[50];
            _itoa_s(expected_nrhs, s1, 10);
            char const s2[] = " input arguments required.";
            size_t const string_lenght = strlen(s1) + 1 + strlen(s2);
            strcat_s(s1, string_lenght, s2);
            mexErrMsgIdAndTxt("MATLAB:Test:nargin", s1);
        }
        else if (nlhs != expected_nlhs)
        {
            char s1[50];
            _itoa_s(expected_nlhs, s1, 10);
            char const s2[] = " output arguments required.";
            size_t const string_lenght = strlen(s1) + 1 + strlen(s2);
            strcat_s(s1, string_lenght, s2);
            mexErrMsgIdAndTxt("MATLAB:Test:nargout", s1);
        }
    }

    int const n_fits = (int)*mxGetPr(prhs[3]);
    int const fit_size = (int)*mxGetPr(prhs[2]);
    int const data_size = n_fits*fit_size;
    int const n_curve_parameters = (int)*mxGetPr(prhs[4]);
    int const functionID = (int)*mxGetPr(prhs[6]);
    real const tolerance = *static_cast<real *>(mxGetData(prhs[7]));

    fcndata_t data_struct;
    data_struct.size = fit_size;
    data_struct.data = new real[data_struct.size];

    real * calculated_function_values = new real[data_struct.size];
    real * jacobian = new real[data_struct.size*n_curve_parameters];
    real * function_parameters;
    int * ipvt = new int[n_curve_parameters];
    int working_array_size = 5 * n_curve_parameters + data_struct.size;
    real * working_array = new real[working_array_size];

    int * info = new int[n_fits];
    int * iterations = new int[n_fits];

    mxArray * mxCurveParameters;
    mxCurveParameters = mxCreateNumericMatrix(1, n_fits*n_curve_parameters, mxSINGLE_CLASS, mxREAL);
    function_parameters = (real*)mxGetData(mxCurveParameters);
    plhs[0] = mxCurveParameters;

    mxArray * mxInfo;
    mxInfo = mxCreateNumericMatrix(1, n_fits, mxINT32_CLASS, mxREAL);
    info = (int*)mxGetData(mxInfo);
    plhs[1] = mxInfo;

    mxArray * mxIterations;
    mxIterations = mxCreateNumericMatrix(1, n_fits, mxINT32_CLASS, mxREAL);
    iterations = (int*)mxGetData(mxIterations);
    plhs[2] = mxIterations;

    for (int ifit = 0; ifit < n_fits; ifit++)
    {
        real const * data_ptr = (real*)mxGetPr(prhs[1]);
        for (int ipixel = 0; ipixel < data_struct.size; ipixel++)
        {
            data_struct.data[ipixel] = real(data_ptr[ifit*data_struct.size + ipixel]);
        }

        real const * x_data_ptr = (real*)mxGetPr(prhs[8]);
        int const x_data_size = x_data_ptr ? data_struct.size : 0;
        data_struct.x_data = new real[x_data_size];
        for (int ipixel = 0; ipixel < x_data_size; ipixel++)
        {
            data_struct.x_data[ipixel] = real(x_data_ptr[ipixel]);
        }

        real const * start_values_ptr = (real*)mxGetPr(prhs[5]);
        real * current_function_parameters = function_parameters + ifit*n_curve_parameters;
        real const * current_start_parameters = (real*)(start_values_ptr)+ifit*n_curve_parameters;
        for (int ipar = 0; ipar < n_curve_parameters; ipar++)
        {
            current_function_parameters[ipar] = current_start_parameters[ipar];
        }

        switch (functionID)
        {
        case GAUSS_1D:
            info[ifit] = __cminpack_func__(lmder1)(
                gauss1d,
                &data_struct,
                data_struct.size,
                n_curve_parameters,
                current_function_parameters,
                calculated_function_values,
                jacobian,
                data_struct.size,
                tolerance,
                ipvt,
                working_array,
                working_array_size);
            break;
        case GAUSS_2D_ELLIPTIC:
            info[ifit] = __cminpack_func__(lmder1)(
                gauss2d_elliptic,
                &data_struct,
                data_struct.size,
                n_curve_parameters,
                current_function_parameters,
                calculated_function_values,
                jacobian,
                data_struct.size,
                tolerance,
                ipvt,
                working_array,
                working_array_size);
            break;
        case GAUSS_2D:
            info[ifit] = __cminpack_func__(lmder1)(
                gauss2d,
                &data_struct,
                data_struct.size,
                n_curve_parameters,
                current_function_parameters,
                calculated_function_values,
                jacobian,
                data_struct.size,
                tolerance,
                ipvt,
                working_array,
                working_array_size);
            break;
        case RAMSEY_FIXED_P:
            info[ifit] = __cminpack_func__(lmder1)(
                calculate_ramsey_fixed_p,
                &data_struct,
                data_struct.size,
                n_curve_parameters,
                current_function_parameters,
                calculated_function_values,
                jacobian,
                data_struct.size,
                tolerance,
                ipvt,
                working_array,
                working_array_size);
            break;
        case RAMSEY_VAR_P:
            info[ifit] = __cminpack_func__(lmder1)(
                calculate_ramsey_var_p,
                &data_struct,
                data_struct.size,
                n_curve_parameters,
                current_function_parameters,
                calculated_function_values,
                jacobian,
                data_struct.size,
                tolerance,
                ipvt,
                working_array,
                working_array_size);
            break;
        case FLETCHER_POWELL:
            info[ifit] = __cminpack_func__(lmder1)(
                calculate_fletcher_powell,
                &data_struct,
                data_struct.size,
                n_curve_parameters,
                current_function_parameters,
                calculated_function_values,
                jacobian,
                data_struct.size,
                tolerance,
                ipvt,
                working_array,
                working_array_size);
            break;
        case BROWN_DENNIS:
            info[ifit] = __cminpack_func__(lmder1)(
                calculate_brown_dennis,
                &data_struct,
                data_struct.size,
                n_curve_parameters,
                current_function_parameters,
                calculated_function_values,
                jacobian,
                data_struct.size,
                tolerance,
                ipvt,
                working_array,
                working_array_size);
            break;
        default:
            break;
        }
        iterations[ifit] = ipvt[0];
    }

    free(data_struct.data);
    free(calculated_function_values);
    free(jacobian);
    free(ipvt);
    free(working_array);

    return;
}