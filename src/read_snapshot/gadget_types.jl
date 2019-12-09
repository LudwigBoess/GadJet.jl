"""
            Types needed for gadget read routines.

    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

mutable struct Header
    npart::Vector{Int32}
    massarr::Vector{Float64}
    time::Float64
    z::Float64
    flag_sfr::Int32
    flag_feedback::Int32
    nall::Vector{UInt32}
    flag_cooling::Int32
    num_files::Int32
    boxsize::Float64
    omega_0::Float64
    omega_l::Float64
    h0::Float64
    flag_stellarage::Int32
    flag_metals::Int32
    npartTotalHighWord::Vector{UInt32}
    flag_entropy_instead_u::Int32
    flag_doubleprecision::Int32
    flag_ic_info::Int32
    lpt_scalingfactor::Float32
    fill::Vector{Int32}

    function Header(npart::Vector{Int32}=Int32.([0,0,0,0,0,0]),
           massarr::Vector{Float64}=zeros(6),
           time::Float64=0.,
           z::Float64=0.,
           flag_sfr::Int32=Int32(0),
           flag_feedback::Int32=Int32(0),
           nall::Vector{UInt32}=UInt32.([0,0,0,0,0,0]),
           flag_cooling::Int32=Int32(0),
           num_files::Int32=Int32(0),
           boxsize::Float64=0.,
           omega_0::Float64=0.,
           omega_l::Float64=0.,
           h0::Float64=0.,
           flag_stellarage::Int32=Int32(0),
           flag_metals::Int32=Int32(0),
           npartTotalHighWord::Vector{UInt32}=UInt32.([0,0,0,0,0,0]),
           flag_entropy_instead_u::Int32=Int32(0),
           flag_doubleprecision::Int32=Int32(0),
           flag_ic_info::Int32=Int32(0),
           lpt_scalingfactor::Float32=Float32(0.),
           fill::Vector{Int32}=Int32.(zeros(12)))

          new(npart,
              massarr,
              time,
              z,
              flag_sfr,
              flag_feedback,
              nall,
              flag_cooling,
              num_files,
              boxsize,
              omega_0,
              omega_l,
              h0,
              flag_stellarage,
              flag_metals,
              npartTotalHighWord,
              flag_entropy_instead_u,
              flag_doubleprecision,
              flag_ic_info,
              lpt_scalingfactor,
              fill)
    end
end

mutable struct Info_Line
    block_name::String
    data_type::DataType
    n_dim::Int32
    is_present::Vector{Int32}

    function Info_Line(block_name="", data_type=Float32, n_dim=Int32(0),
                        is_present=Int32.(zeros(6)))

        new(block_name, data_type, n_dim, is_present)
    end
end

mutable struct part_single
      ptype::Int64
      pos::Vector{Float32}
      vel::Vector{Float32}
      id::UInt32
      m::Float32
      U::Float32
      rho::Float32
      sml::Float32

      function part_single(ptype=-1, pos=Float32.(zeros(3)), vel=Float32.(zeros(3)),
                        id=UInt32(0), m=Float32(0.), U=Float32(0.), rho=Float32(0.),
                        sml=Float32(0.))

             new(ptype, pos, vel, id, m, U, rho, sml)

      end
end

mutable struct part_double
      ptype::Int64
      pos::Vector{Float64}
      vel::Vector{Float64}
      id::UInt32
      m::Float64
      U::Float64
      rho::Float64
      sml::Float64

      function part_double(ptype=-1, pos=zeros(3), vel=zeros(3),
                        id=UInt32(0), m=0., U=0., rho=0.,
                        sml=0.)

             new(ptype, pos, vel, id, m, U, rho, sml)

      end
end
