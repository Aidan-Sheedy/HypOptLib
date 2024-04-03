import signac

project = signac.init_project()

for temperature in [0]:
    sp = {"T":          temperature,
          "numx":       64,
          "numy":       32,
          "numz":       32,
          "dt":         0.1,
          "maxitr":     5000,#00,
          "saveFreq":   1000,
          "filePath":   "hypopt_output_t"+str(temperature)+".h5"
          }
    print(sp)
    job = project.open_job(sp).init()
