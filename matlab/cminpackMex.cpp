#include <iostream>
#include "mex.h"
#include <stdlib.h>
#include <vector>

#define real __cminpack_real__
#include "..\cminpack.h"

#define GAUSS_1D 0
#define GAUSS_2D 1
#define GAUSS_2D_ELLIPTIC 2
#define GAUSS_2D_ROTATED 3
#define CAUCHY_2D_ELLIPTIC 4

typedef struct  {
    int size;
    real * data;
} fcndata_t;

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

void mexFunction(
    int          nlhs,
    mxArray      *plhs[],
    int          nrhs,
    mxArray const *prhs[])
{
    int expected_nrhs = 8;
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
    float const tolerance = (float)*mxGetPr(prhs[7]);

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
    mxCurveParameters = mxCreateNumericMatrix(1, n_fits*n_curve_parameters, mxDOUBLE_CLASS, mxREAL);
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
        double const * data_ptr = mxGetPr(prhs[1]);
        for (int ipixel = 0; ipixel < data_struct.size; ipixel++)
        {
            data_struct.data[ipixel] = real(data_ptr[ifit*data_struct.size + ipixel]);
        }

        double const * start_values_ptr = (mxGetPr(prhs[5]));
        real * current_function_parameters = function_parameters + ifit*n_curve_parameters;
        real const * current_start_parameters = (real*)(start_values_ptr)+ifit*n_curve_parameters;
        for (int ipar = 0; ipar < n_curve_parameters; ipar++)
        {
            current_function_parameters[ipar] = current_start_parameters[ipar];
        }

        switch (functionID)
        {
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