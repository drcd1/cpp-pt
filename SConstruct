env = Environment(CPPPATH = ['src', 'third-party'])
env.Append(CCFLAGS = '/O2 /openmp /MT /EHsc')
#env.Append(LINKFLAGS = '/DEBUG')
env.Program("build/renderer","src/main.cpp")