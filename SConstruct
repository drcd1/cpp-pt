env = Environment(CPPPATH = ['src', 'third-party'])

debug = ARGUMENTS.get('debug', 0)

if (int(debug)):
    env.Append(CCFLAGS = '/Od /openmp /MT /EHsc')
    env.Append(PDB = 'build/main.pdb')
    env.Append(CCPDBFLAGS = '/Zi')
    env.Append(LINKFLAGS = '/DEBUG')
else:
    env.Append(CCFLAGS = "/O2 /openmp /MT /EHsc")

env.Program("build/renderer","src/main.cpp")

