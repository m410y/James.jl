function peak_hkl(experiment::Experiment, peak::Peak)
    setting = peak.frame.setting
    coord = peak.coord
    peak_model = model(experiment, setting)
    round.(reflex(peak_model, coord), RoundNearest)
end

function experiment_correction(experiment::Experiment, p)
    detector_trans = AffineMap(RotationVec(p[1:3]...), p[4:6]u"mm")
    detector0 = experiment.detector.object
    detector = detector_trans(detector0)
    crystal_trans = AffineMap(RotationVec(p[7:9]...), p[10:12]u"mm")
    crystal0 = experiment.sample.object
    crystal = crystal_trans(crystal0)
    Experiment(
        experiment.spectrum,
        Motorized(experiment.detector.goniometer, detector),
        Motorized(experiment.sample.goniometer, crystal),
    )
end

function refine_experiment(experiment::Experiment, peaks::AbstractArray{Peak})
    coords = [peak.coord for peak in peaks]
    settings = [peak.frame.setting for peak in peaks]
    hkls = [peak_hkl(experiment, peak) for peak in peaks]
    func =
        p -> begin
            corrected_experiment = experiment_correction(experiment, p)
            [
                coord(model(corrected_experiment, setting), hkl) for
                (setting, hkl) in zip(settings, hkls)
            ]
        end
    p0 = zeros(12)
    pmin = fill(-0.4, 12)
    pmax = fill(0.4, 12)
    lsq = OptimizationFunction(
        (p, _) -> sum(LinearAlgebra.norm2, coords .- func(p)),
        AutoFiniteDiff(),
    )
    fit = solve(OptimizationProblem(lsq, p0, lb = pmin, ub = pmax), Optimization.LBFGS())
    @show fit
    experiment_correction(experiment, fit)
end
