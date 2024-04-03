from flow import FlowProject
import cantilevered_beam


class Project(FlowProject):
    pass

class HypOptCAC(environment.DefaultSlurmSnevironment):
    hostname_pattern=r"cac.*"
    template = "HypOptCAC.sh"

def is_completed(job):
    return job.isfile("Finished.txt")

@Project.post(is_completed)
# @Project.operation
@FlowProject.operation(directives={"nranks": 20})
def simulate_hypoptlib(job):
    with job:
        beam = cantilevered_beam.CantileveredBeam(
            [job.sp.numx, job.sp.numy, job.sp.numz],
            job.sp.dt,
            job.sp.maxitr,
            job.sp.T,
            job.sp.saveFreq,
            job.sp.filePath
        )
        beam.startRun([0,job.sp.maxitr])
        with open(job.fn("Finished.txt"), "w") as file:
            file.write("finished\n")

if __name__ == "__main__":
    Project().main()
