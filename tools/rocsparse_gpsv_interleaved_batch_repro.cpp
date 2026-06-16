#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>

#if __has_include(<hipsparse/hipsparse.h>)
#include <hipsparse/hipsparse.h>
#elif __has_include(<hipsparse.h>)
#include <hipsparse.h>
#else
#error "Could not find hipSPARSE headers"
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
constexpr int block_size = 256;

[[noreturn]] void die(const std::string& message)
{
    std::cerr << "error: " << message << "\n";
    std::exit(1);
}

void check_hip(hipError_t status, const char* where)
{
    if(status != hipSuccess)
    {
        std::cerr << "HIP error in " << where << ": " << hipGetErrorString(status) << " ("
                  << static_cast<int>(status) << ")\n";
        std::exit(2);
    }
}

void check_sparse(hipsparseStatus_t status, const char* where)
{
    if(status != HIPSPARSE_STATUS_SUCCESS)
    {
        std::cerr << "hipSPARSE error in " << where << ": status=" << static_cast<int>(status)
                  << "\n";
        std::exit(3);
    }
}

std::uint64_t parse_u64(const char* text, const char* name)
{
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if(end == text || *end != '\0')
    {
        die(std::string("invalid integer for ") + name + ": " + text);
    }
    return static_cast<std::uint64_t>(value);
}

__host__ __device__ hipDoubleComplex zmake(double real, double imag = 0.0)
{
    return make_hipDoubleComplex(real, imag);
}

__global__ void init_penta_systems(hipDoubleComplex* ds,
                                   hipDoubleComplex* dl,
                                   hipDoubleComplex* d,
                                   hipDoubleComplex* du,
                                   hipDoubleComplex* dw,
                                   hipDoubleComplex* x,
                                   int               m,
                                   int               batch_count,
                                   std::size_t       total)
{
    const std::size_t tid = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if(tid >= total)
    {
        return;
    }

    const int row = static_cast<int>(tid / static_cast<std::size_t>(batch_count));
    const int b   = static_cast<int>(tid - static_cast<std::size_t>(row) * batch_count);

    const double tweak = static_cast<double>(b & 15) * 1.0e-8;
    const double a2    = (row >= 2) ? 0.005 : 0.0;
    const double a1    = (row >= 1) ? 0.020 : 0.0;
    const double a0    = 2.000 + tweak;
    const double c1    = (row + 1 < m) ? 0.015 : 0.0;
    const double c2    = (row + 2 < m) ? 0.004 : 0.0;

    ds[tid] = zmake(a2);
    dl[tid] = zmake(a1);
    d[tid]  = zmake(a0);
    du[tid] = zmake(c1);
    dw[tid] = zmake(c2);

    // Exact solution is x_true = 1 + 0i for every system.
    x[tid] = zmake(a2 + a1 + a0 + c1 + c2);
}

__global__ void fill_complex(hipDoubleComplex* x, std::size_t total, hipDoubleComplex value)
{
    const std::size_t tid = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if(tid < total)
    {
        x[tid] = value;
    }
}

__global__ void validate_solution(const hipDoubleComplex* x,
                                  std::size_t             total,
                                  double*                 block_max_error,
                                  unsigned long long*     block_bad_count)
{
    __shared__ double             smax[block_size];
    __shared__ unsigned long long sbad[block_size];

    const int         lane = threadIdx.x;
    const std::size_t tid  = static_cast<std::size_t>(blockIdx.x) * blockDim.x + lane;

    double             err = 0.0;
    unsigned long long bad = 0;

    if(tid < total)
    {
        const hipDoubleComplex value = x[tid];
        const double           real  = value.x;
        const double           imag  = value.y;
        if(!isfinite(real) || !isfinite(imag))
        {
            bad = 1;
            err = INFINITY;
        }
        else
        {
            err = hypot(real - 1.0, imag);
        }
    }

    smax[lane] = err;
    sbad[lane] = bad;
    __syncthreads();

    for(int stride = blockDim.x / 2; stride > 0; stride >>= 1)
    {
        if(lane < stride)
        {
            smax[lane] = fmax(smax[lane], smax[lane + stride]);
            sbad[lane] += sbad[lane + stride];
        }
        __syncthreads();
    }

    if(lane == 0)
    {
        block_max_error[blockIdx.x] = smax[0];
        block_bad_count[blockIdx.x] = sbad[0];
    }
}

__global__ void validate_zero(const hipDoubleComplex* x,
                              std::size_t             total,
                              double*                 block_max_abs,
                              unsigned long long*     block_bad_count)
{
    __shared__ double             smax[block_size];
    __shared__ unsigned long long sbad[block_size];

    const int         lane = threadIdx.x;
    const std::size_t tid  = static_cast<std::size_t>(blockIdx.x) * blockDim.x + lane;

    double             err = 0.0;
    unsigned long long bad = 0;

    if(tid < total)
    {
        const hipDoubleComplex value = x[tid];
        const double           real  = value.x;
        const double           imag  = value.y;
        if(!isfinite(real) || !isfinite(imag))
        {
            bad = 1;
            err = INFINITY;
        }
        else
        {
            err = hypot(real, imag);
        }
    }

    smax[lane] = err;
    sbad[lane] = bad;
    __syncthreads();

    for(int stride = blockDim.x / 2; stride > 0; stride >>= 1)
    {
        if(lane < stride)
        {
            smax[lane] = fmax(smax[lane], smax[lane + stride]);
            sbad[lane] += sbad[lane + stride];
        }
        __syncthreads();
    }

    if(lane == 0)
    {
        block_max_abs[blockIdx.x]  = smax[0];
        block_bad_count[blockIdx.x] = sbad[0];
    }
}

struct Validation
{
    double             max_value = 0.0;
    unsigned long long bad_count = 0;
};

Validation collect_validation(double*             d_max,
                              unsigned long long* d_bad,
                              std::size_t         blocks64,
                              const char*         label)
{
    check_hip(hipGetLastError(), label);
    check_hip(hipDeviceSynchronize(), label);

    std::vector<double>             h_max(blocks64);
    std::vector<unsigned long long> h_bad(blocks64);
    check_hip(hipMemcpy(h_max.data(), d_max, h_max.size() * sizeof(double), hipMemcpyDeviceToHost),
              "copy validation max");
    check_hip(hipMemcpy(h_bad.data(), d_bad, h_bad.size() * sizeof(unsigned long long), hipMemcpyDeviceToHost),
              "copy validation bad");

    Validation result;
    for(std::size_t i = 0; i < blocks64; ++i)
    {
        result.max_value = std::max(result.max_value, h_max[i]);
        result.bad_count += h_bad[i];
    }
    return result;
}

Validation run_validate_zero(const hipDoubleComplex* x, std::size_t total)
{
    const std::size_t blocks64 = (total + block_size - 1) / block_size;
    if(blocks64 > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        die("validation grid is too large for a one-dimensional launch");
    }
    const int blocks = static_cast<int>(blocks64);

    double*             d_max = nullptr;
    unsigned long long* d_bad = nullptr;
    check_hip(hipMalloc(&d_max, blocks64 * sizeof(double)), "hipMalloc validation max");
    check_hip(hipMalloc(&d_bad, blocks64 * sizeof(unsigned long long)), "hipMalloc validation bad");

    validate_zero<<<blocks, block_size>>>(x, total, d_max, d_bad);
    Validation result = collect_validation(d_max, d_bad, blocks64, "validate_zero");

    check_hip(hipFree(d_max), "hipFree validation max");
    check_hip(hipFree(d_bad), "hipFree validation bad");
    return result;
}

Validation run_validate_solution(const hipDoubleComplex* x, std::size_t total)
{
    const std::size_t blocks64 = (total + block_size - 1) / block_size;
    if(blocks64 > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        die("validation grid is too large for a one-dimensional launch");
    }
    const int blocks = static_cast<int>(blocks64);

    double*             d_max = nullptr;
    unsigned long long* d_bad = nullptr;
    check_hip(hipMalloc(&d_max, blocks64 * sizeof(double)), "hipMalloc validation max");
    check_hip(hipMalloc(&d_bad, blocks64 * sizeof(unsigned long long)), "hipMalloc validation bad");

    validate_solution<<<blocks, block_size>>>(x, total, d_max, d_bad);
    Validation result = collect_validation(d_max, d_bad, blocks64, "validate_solution");

    check_hip(hipFree(d_max), "hipFree validation max");
    check_hip(hipFree(d_bad), "hipFree validation bad");
    return result;
}

void print_bytes(const char* label, std::size_t bytes)
{
    const double gib = static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
    std::cout << label << bytes << " bytes (" << std::fixed << std::setprecision(3) << gib
              << " GiB)\n";
}

void usage(const char* argv0)
{
    std::cerr << "usage: " << argv0 << " [--m N] [--batch-count N] [--repeat N]\n"
              << "       " << argv0 << " --memset-only [--m N] [--batch-count N]\n"
              << "\n"
              << "Defaults match the nx=511, nz=257, ny=520-ish endpoint-Schur batch:\n"
              << "  --m 258 --batch-count 527872\n"
              << "Try --m 253 for the just-below case.\n";
}

} // namespace

int main(int argc, char** argv)
{
    int  m                  = 258;
    int  batch_count        = 527872;
    int  repeat             = 1;
    bool run_memset_check   = true;
    bool run_sparse_check   = true;

    for(int i = 1; i < argc; ++i)
    {
        if(std::strcmp(argv[i], "--m") == 0 && i + 1 < argc)
        {
            m = static_cast<int>(parse_u64(argv[++i], "--m"));
        }
        else if((std::strcmp(argv[i], "--batch-count") == 0 || std::strcmp(argv[i], "--batch") == 0)
                && i + 1 < argc)
        {
            batch_count = static_cast<int>(parse_u64(argv[++i], "--batch-count"));
        }
        else if(std::strcmp(argv[i], "--repeat") == 0 && i + 1 < argc)
        {
            repeat = static_cast<int>(parse_u64(argv[++i], "--repeat"));
        }
        else if(std::strcmp(argv[i], "--skip-memset") == 0)
        {
            run_memset_check = false;
        }
        else if(std::strcmp(argv[i], "--memset-only") == 0)
        {
            run_sparse_check = false;
        }
        else if(std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0)
        {
            usage(argv[0]);
            return 0;
        }
        else
        {
            usage(argv[0]);
            return 1;
        }
    }

    if(m < 3 || batch_count < 1 || repeat < 1)
    {
        die("m must be >= 3, batch-count >= 1, repeat >= 1");
    }

    const std::size_t total = static_cast<std::size_t>(m) * static_cast<std::size_t>(batch_count);
    if(total > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        std::cout << "warning: m*batch_count exceeds INT_MAX; this intentionally goes beyond "
                     "rocSPARSE's 32-bit element indexing\n";
    }

    int device = 0;
    check_hip(hipGetDevice(&device), "hipGetDevice");
    hipDeviceProp_t props{};
    check_hip(hipGetDeviceProperties(&props, device), "hipGetDeviceProperties");

    std::size_t free_bytes = 0;
    std::size_t total_bytes = 0;
    check_hip(hipMemGetInfo(&free_bytes, &total_bytes), "hipMemGetInfo");

    std::cout << "device=" << device << " name=\"" << props.name << "\"\n";
    std::cout << "m=" << m << " batch_count=" << batch_count << " total_complex=" << total
              << "\n";
    print_bytes("one complex array: ", total * sizeof(hipDoubleComplex));
    print_bytes("six matrix/rhs arrays: ", 6 * total * sizeof(hipDoubleComplex));
    print_bytes("free device memory: ", free_bytes);

    const std::size_t blocks64 = (total + block_size - 1) / block_size;
    if(blocks64 > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        die("initialization grid is too large for a one-dimensional launch");
    }
    const dim3 blocks(static_cast<unsigned int>(blocks64));
    const dim3 threads(block_size);

    if(run_memset_check)
    {
        std::cout << "\nraw hipMemsetAsync check\n";
        hipDoubleComplex* scratch = nullptr;
        check_hip(hipMalloc(&scratch, total * sizeof(hipDoubleComplex)), "hipMalloc memset scratch");
        fill_complex<<<blocks, threads>>>(scratch, total, zmake(3.0, -4.0));
        check_hip(hipGetLastError(), "fill_complex launch");
        check_hip(hipDeviceSynchronize(), "fill_complex");

        check_hip(hipMemsetAsync(scratch, 0, total * sizeof(hipDoubleComplex), nullptr),
                  "hipMemsetAsync scratch");
        check_hip(hipDeviceSynchronize(), "hipMemsetAsync scratch sync");

        const Validation zero = run_validate_zero(scratch, total);
        std::cout << "memset max_abs=" << std::scientific << zero.max_value
                  << " bad_count=" << zero.bad_count << "\n";
        check_hip(hipFree(scratch), "hipFree memset scratch");
    }

    if(!run_sparse_check)
    {
        return 0;
    }

    hipsparseHandle_t handle = nullptr;
    check_sparse(hipsparseCreate(&handle), "hipsparseCreate");

    hipDoubleComplex *ds = nullptr, *dl = nullptr, *d = nullptr, *du = nullptr, *dw = nullptr, *x = nullptr;
    check_hip(hipMalloc(&ds, total * sizeof(hipDoubleComplex)), "hipMalloc ds");
    check_hip(hipMalloc(&dl, total * sizeof(hipDoubleComplex)), "hipMalloc dl");
    check_hip(hipMalloc(&d, total * sizeof(hipDoubleComplex)), "hipMalloc d");
    check_hip(hipMalloc(&du, total * sizeof(hipDoubleComplex)), "hipMalloc du");
    check_hip(hipMalloc(&dw, total * sizeof(hipDoubleComplex)), "hipMalloc dw");
    check_hip(hipMalloc(&x, total * sizeof(hipDoubleComplex)), "hipMalloc x");

    init_penta_systems<<<blocks, threads>>>(ds, dl, d, du, dw, x, m, batch_count, total);
    check_hip(hipGetLastError(), "init_penta_systems launch");
    check_hip(hipDeviceSynchronize(), "init_penta_systems");

    std::size_t buffer_size = 0;
    check_sparse(hipsparseZgpsvInterleavedBatch_bufferSizeExt(
                     handle, 0, m, ds, dl, d, du, dw, x, batch_count, &buffer_size),
                 "hipsparseZgpsvInterleavedBatch_bufferSizeExt");
    print_bytes("hipSPARSE temp buffer: ", buffer_size);

    void* buffer = nullptr;
    check_hip(hipMalloc(&buffer, buffer_size), "hipMalloc gpsv buffer");

    for(int iter = 0; iter < repeat; ++iter)
    {
        if(iter > 0)
        {
            init_penta_systems<<<blocks, threads>>>(ds, dl, d, du, dw, x, m, batch_count, total);
            check_hip(hipGetLastError(), "reinit_penta_systems launch");
            check_hip(hipDeviceSynchronize(), "reinit_penta_systems");
        }

        std::cout << "\nhipSPARSE solve iteration " << (iter + 1) << "\n";
        check_sparse(hipsparseZgpsvInterleavedBatch(
                         handle, 0, m, ds, dl, d, du, dw, x, batch_count, buffer),
                     "hipsparseZgpsvInterleavedBatch");
        check_hip(hipDeviceSynchronize(), "hipsparseZgpsvInterleavedBatch sync");

        const Validation solution = run_validate_solution(x, total);
        std::cout << "solution max_error=" << std::scientific << solution.max_value
                  << " bad_count=" << solution.bad_count << "\n";
    }

    check_hip(hipFree(buffer), "hipFree gpsv buffer");
    check_hip(hipFree(x), "hipFree x");
    check_hip(hipFree(dw), "hipFree dw");
    check_hip(hipFree(du), "hipFree du");
    check_hip(hipFree(d), "hipFree d");
    check_hip(hipFree(dl), "hipFree dl");
    check_hip(hipFree(ds), "hipFree ds");
    check_sparse(hipsparseDestroy(handle), "hipsparseDestroy");

    return 0;
}
