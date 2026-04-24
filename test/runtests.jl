using Test
using SplitVOF
using StaticArrays
using VOFTools
using LinearAlgebra

include("test_init.jl")
include("test_reconstruction.jl")
include("test_translation.jl")
include("test_rotation.jl")
include("test_lie_splitting.jl")
include("test_3d_init.jl")
include("test_3d_reconstruction.jl")
include("test_3d_advection.jl")
