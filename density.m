function rho = density(n_e, A, Z)
    % calculates mass density from electron density (e- per Angstrom^3),
    % molar mass (g/mol), and electrons per mole

    N_A = 6.0221367e23; % mol^-1, Avogadro constant

    rho = (n_e .* A) ./ (N_A .* Z) .* 1e24; % g/cm^3
end