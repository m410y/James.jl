const SpaceUnit = u"mm"
const ReciprocalUnit = u"eV"

reciprocal_convert(unit, E) =
    unit != ReciprocalUnit ? ustrip(uconvert(unit, E * ReciprocalUnit, Spectral())) : E

space_convert(unit, L) = unit != SpaceUnit ? ustrip(uconvert(unit, L * SpaceUnit)) : L

angle_norm(angle::Real) = abs(rem2pi(angle, RoundNearest))
