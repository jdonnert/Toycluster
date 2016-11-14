module Toycluster  # experimental

const Msol = 1.989e33

type Cluster

	beta::Float64
	rc::Float64
	rho0_gas::Float64
	rho0_dm::Float64
	M200::Float64
	R200::Float64
	rs::Float64
	c_NFW::Float64
	rcut::Float64
	bf::Float64
	z::Float64

	function Cluster(M200::Float64, beta::Float64; rc=-1, c_NFW=-1, rcut=-1, bf=0.17, z=0, coolCore=false)

		r200 = R200(M200*1e10)
	
		if c_NFW == -1
			c_nfw = c_NFW(M200*1e10, z)
		end

		rs = r200/c_nfw
		
		if rc == -1
			if coolCore == true
				rc = rs/9
			else
				rc = rs/3
			end
		end

		if rcut == -1
			rcut = r200
		end

		cluster = new(beta, rc, 1, 1, M200, r200, rs, c_NFW, rcut)

		return cluster
	end

end

function R200(mass; delta=103.69, rho_crit=1.10744e-29)

	return  (mass*Msol ./ (delta*rho_crit * 4*pi/3))^(1/3) / kpc2cm
end

function c_NFW(mass, z)

	A = 5.74 # Duffy+ 2008
	B = -0.097
	C = -0.47
	Mpivot = 2e12 / 0.7

	return A .* (mass./Mpivot).^(B) .* (1+z)^C
	
end

function Mgas(gc::Cluster)

	return gc.M200 * gc.bf
end

function DMDensityProfile(gc::Cluster, r)
	
	rr = r./gc.rs

	return gc.rho0_dm ./ (rr .* (1.+rr).^2)
end

function GasDensityProfile(gc::Cluster, r)

	rho_beta = gc.rho0_gas .* (1 .+ r.^2 ./ gc.rc^2 ).^(-3/2*gc.beta)

	return rho_beta./(1 .+ r.^3 ./ gc.rcut.^3)
end

function DMMassProfile(gc::Cluster, r)

	return 4*pi * gc.rho0_dm * gc.rs^3 .* (log( (r.+ gc.rs)/gc.rs) - r./(gc.rs.+r))
end

function MProf_Integrant(gc::Cluster, r)

	return 4 * pi .* r.^2 .* GasDensityProfile(r, gc)
end

function GasMassProfile(gc::Cluster, r)
		
	nStep = 256

	N = max(length(r),1)

	Mr = zeros(N)

	for i = 1:N

		x = logspace(-3, log10(r[i]), nStep)

		dx = x[2:nStep] - x[1:nStep-1]

		f = MProf_Integrant(x, gc)

		Mr[i] = sum(dx .* 0.5 .* (f[2:nStep] + f[1:nStep-1]))
	end
	
	return Mr
end

function MagneticFieldProfile(gc::Cluster, r)

	return 1e-5 * (GasDensityProfile(r,gc)/gc.rho0_gas).^(0.5)
end

function UIntegrant(gc::Cluster, r)

	Mr = GasMassProfile(r, gc) .+ DMMassProfile(r, gc)

	return GasDensityProfile(r, gc) .* Mr ./ r.^2

end

function InternalEnergy(gc::Cluster, r::Array)

	nStep = 512

	N = max(length(r),1)

	U = Array(Float64, N)

	for i = 1:N
		
		x = logspace(log10(r[i]), 10, nStep)
		
		dx = x[2:nStep] - x[1:nStep-1]

		f = UIntegrant(x, gc)

		U[i] = sum(dx .* 0.5 .* (f[2:nStep] + f[1:nStep-1]))
	end

	rho_Gas = GasDensityProfile(r, gc)

	U .*=  43007./( (gamma-1) .* rho_Gas )

	return U
end

function SoundSpeed(gc::Cluster, r::Array)

	nStep = 512

	N = max(length(r),1)

	U = InternalEnergy(r, gc)

	return sqrt( U *gamma * (gamma-1) )
end



end #module
