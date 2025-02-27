function peak_hkl(experiment::Experiment, peak::Peak)
    setting = peak.frame.setting
    coord = peak.coord
    model = model(experiment, setting)
    round.(reflex(model, coord), RoundNearest)
end

function experiment_correction(experiment::Experiment, p)
    crystal = SingleCrystal(
        Vec3(p[1:3]),
        Mat3(p[4:12])
    )
    Experiment(
        experiment.spectrum,
        experiment.detector,
        Motorized(experiment.sample.goniometer, crystal)
    )
end

function refine_experiment(experiment::Experiment, peaks::AbstractArray{Peak})
    coords = [peak.coord for peak in peaks]
    settings = [peak.frame.setting for peak in peaks]
    hkls = [peak_hkl(experiment, peak) for peak in peaks]
    func = p -> begin
        corrected_experiment = experiment_correction(experiment, p)
        [coord(model(corrected_experiment, setting), hkl) for (setting, hkl) in zip(settings, hkls)]
    end
    lsq = OptimizationFunction((p, _) -> sum(LinearAlgebra.norm2, coords .- func(p)), AutoZygote())
    fit = solve(OptimizationProblem(lsq, experiment), Optimization.LBFGS())
    experiment_correction(experiment, fit)
end
