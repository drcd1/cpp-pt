env = Environment(CPPPATH = ['src', 'third-party'])
env.Append(CCFLAGS = '/Od /openmp /MT /EHsc')
env.Append(PDB = 'build/main.pdb')
env.Append(CCPDBFLAGS = '/Zi')
env.Append(LINKFLAGS = '/DEBUG')
env.Program("build/renderer","src/main.cpp")

