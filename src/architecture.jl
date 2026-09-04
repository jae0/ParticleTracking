"""
    architecture.jl

Computational architecture abstraction and hardware resolution (CPU vs NVIDIA GPU/CUDA)
for Oceananigans.jl hydrodynamic simulations and Lagrangian particle tracking.
"""

using Oceananigans
using Oceananigans.Architectures: AbstractArchitecture, CPU, GPU

"""
    resolve_architecture(
        arch::Union{Symbol, String, Bool, AbstractArchitecture} = :cpu;
        fallback_to_cpu::Bool = false
    ) -> AbstractArchitecture

Resolve and instantiate the computational architecture (`CPU()` or `GPU(...)`) for
Oceananigans hydrodynamic simulations and particle tracking routines.

# Mathematical & System Architecture
Oceananigans supports multi-threaded CPU execution and massively parallel NVIDIA CUDA GPU
acceleration through Julia's GPUCompiler / CUDA.jl backend. On GPUs, array allocations
use unified device memory and compute kernels run asynchronously across streaming multiprocessors.

```math
\\text{Architecture} \\in \\{ \\text{CPU}(), \\text{GPU}(\\text{device}) \\}
```

# Arguments
- `arch`: Architecture descriptor. Supported values:
  - `:cpu`, `"cpu"`, or `false`: Explicit multi-threaded CPU architecture (`CPU()`).
  - `:gpu`, `:cuda`, `"gpu"`, `"cuda"`, or `true`: NVIDIA CUDA GPU architecture (`GPU()`).
  - `Oceananigans.Architectures.AbstractArchitecture`: Pre-constructed architecture instance.
- `fallback_to_cpu::Bool`: If `true`, issues an informative warning and returns `CPU()`
  when GPU is requested on hardware without functional CUDA drivers. If `false` (default),
  raises an informative `ErrorException`.

# Returns
- `AbstractArchitecture`: Concrete instantiated architecture (`CPU()` or `GPU(...)`).

# Throws
- `ErrorException`: When GPU is requested without functional CUDA runtime and `fallback_to_cpu=false`.
"""
function resolve_architecture(
    arch::Union{Symbol, String, Bool, AbstractArchitecture} = :cpu;
    fallback_to_cpu::Bool = false
)::AbstractArchitecture
    if arch isa AbstractArchitecture
        return arch
    end

    is_gpu_req = (arch === true) || (arch === :gpu) || (arch === :cuda) ||
                 (arch isa String && lowercase(arch) in ["gpu", "cuda"])

    if !is_gpu_req
        return CPU()
    end

    # Check CUDA availability
    cuda_functional = false
    if !isdefined(Main, :CUDA)
        try
            # Attempt to dynamically import CUDA if installed in the environment
            if Base.find_package("CUDA") !== nothing
                @eval Main import CUDA
            end
        catch err
            @debug "Dynamic CUDA import skipped: $(err)"
        end
    end

    try
        if isdefined(Main, :CUDA)
            cuda_functional = Main.CUDA.functional()
        elseif isdefined(Oceananigans, :CUDA)
            cuda_functional = Oceananigans.CUDA.functional()
        end
    catch
        cuda_functional = false
    end

    if cuda_functional
        try
            if isdefined(Oceananigans, :GPU)
                return Oceananigans.GPU()
            elseif isdefined(Oceananigans.Architectures, :GPU)
                return Oceananigans.Architectures.GPU()
            end
        catch err
            @warn "Failed to construct Oceananigans.GPU device: $(err)"
        end
    end

    if fallback_to_cpu
        @warn "CUDA GPU hardware was requested, but no functional NVIDIA CUDA environment " *
              "was detected (ensure the CUDA.jl package is installed and NVIDIA drivers " *
              "are accessible). Falling back to CPU()."
        return CPU()
    else
        error(
            "CUDA GPU acceleration was requested (`arch = $(arch)`), but no functional " *
            "NVIDIA CUDA driver or GPU device was detected on this system.\n" *
            "To execute on CPU, set `architecture = :cpu` (or `--cpu` CLI flag).\n" *
            "To enable automatic fallback, set `fallback_to_cpu = true`."
        )
    end
end
