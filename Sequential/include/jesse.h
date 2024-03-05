void calc_force() {
  potential_sum = 0.0;

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {

      force_x = 0.0;
      force_y = 0.0;
      force_z = 0.0;

      vector<double> distance = particle_distance(positions[i], positions[j]);
      vector<double> micdistance = minimum_image_condition(distance);
      force_x += 48.0 * pow(magnitude(micdistance), -8.0) *
                 (pow(magnitude(micdistance), -6.0) - 0.5) * micdistance.at(0);
      force_y += 48.0 * pow(magnitude(micdistance), -8.0) *
                 (pow(magnitude(micdistance), -6.0) - 0.5) * micdistance.at(1);
      force_z += 48.0 * pow(magnitude(micdistance), -8.0) *
                 (pow(magnitude(micdistance), -6.0) - 0.5) * micdistance.at(2);
      // A mathematician told me this code is bad and I should feel bad
      vector<double> distance = particle_distance(positions[i], positions[j]);
      vector<double> micdistance = minimum_image_condition(distance);
      potential_sum += potential_energy(magnitude(micdistance));
    }
  }
}
