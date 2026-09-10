
using Oceananigans
using Oceananigans.Architectures: AbstractArchitecture, CPU, GPU
using CUDA

"""
    resolve_architecture(arch::Union{Symbol, String, Bool, AbstractArchitecture} = :cpu; 
                         fallback_to_cpu::Bool = false) -> AbstractArchitecture

Resolve and instantiate the computational architecture (`CPU()` or `GPU(...)`) for 
Oceananigans hydrodynamic simulations and particle tracking routines, with world-age safety for Julia 1.12+.
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

    # Check CUDA availability safely avoiding world-age issues under Julia 1.12+ and Revise
    cuda_functional = false
    try
        if isdefined(Main, :CUDA) && isdefined(Main.CUDA, :functional)
            cuda_functional = Base.invokelatest(Main.CUDA.functional)
        elseif isdefined(Oceananigans, :CUDA) && isdefined(Oceananigans.CUDA, :functional)
            cuda_functional = Base.invokelatest(Oceananigans.CUDA.functional)
        else
            # Attempt dynamic check if CUDA package is loaded or loadable
            cuda_functional = Base.invokelatest(CUDA.functional)
        end
    catch err
        @debug "CUDA availability check encountered an error: $(err)"
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