
set(SimulatorObjects
  "Settings.cpp"
  "SPHSolver.cpp"
  "Kernels.cpp"
  "HashTable.cpp"
  "Reorderer.cpp"
)

add_library (Simulator ${SimulatorObjects})
