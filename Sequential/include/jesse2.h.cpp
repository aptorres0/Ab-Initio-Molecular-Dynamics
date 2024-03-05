void calc_force() {
  potential_sum = 0.0;

  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      vector<double> distance = particle_distance(positions[i], positions[j]);
      vector<double> micdistance = minimum_image_condition(distance);
      vector<double> distance = particle_distance(positions[i], positions[j]);
      vector<double> micdistance = minimum_image_condition(distance);
      potential_sum += potential_energy(magnitude(micdistance));
    }
  }
}