#ifndef ACTUAL_EOS_H
#define ACTUAL_EOS_H

#include <string>
#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <cmath>

using namespace eos_rp;

// This is a constant gamma equation of state, using an ideal gas.
//
// The gas may either be completely ionized or completely neutral.
//
// The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
// expression for entropy is only valid for an ideal MONATOMIC gas
// (gamma = 5/3).

const std::string eos_name = "gamma_law_2T";

inline
void actual_eos_init() {

  // make sure we have the needed aux data

  // ensure that there are only 2 species, one with Z = 0

  // ensure that the species have the same atomic mass

  // constant ratio of specific heats
  if (eos_gamma <= 0.0) {
    amrex::Error("gamma_const cannot be < 0");
  }

}


template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_input_valid (I input)
{
  static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

  bool valid = true;

  if (input == eos_input_th) {
      valid = false;
  } else if (input == eos_input_rh) {
      valid = false;
  } else if (input == eos_input_tp) {
      valid = false;
  } else if (input == eos_input_ps) {
      valid = false;
  } else if (input == eos_input_ph) {
      valid = false;
  } else if (input == eos_input_th) {
      valid = false;
  }

  return valid;
}


template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_eos (I input, T& state)
{
    static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

    // Get the mass of a nucleon from m_u.
    const amrex::Real m_nucleon = C::m_u;

    if constexpr (has_xn<T>::value) {
        amrex::Real sum = 0.0;
        for (int n = 0; n < NumSpec; n++) {
            sum += (zion[n] + 1.0) * state.xn[n] * aion_inv[n];
        }
        state.mu = 1.0 / sum;
    }

    // Find the fraction of electrons
    amrex::Real x_e = 0.0;
    for (int n = 0; n < NumSpec; n++) {
        x_e += zion[n] * state.xn[n] * aion_inv[n];
    }

    switch (input)
    {

      case eos_input_rt:

          // dens, temp, xmass are inputs

          // We don't need to do anything here
          break;

      case eos_input_rh:

        // dens, enthalpy, and xmass are inputs

#ifndef AMREX_USE_GPU
          amrex::Error("EOS: eos_input_tp has not been implemented gamma_law_2T EOS.");
#endif

          break;

      case eos_input_tp:

          // temp, pres, and xmass are inputs

#ifndef AMREX_USE_GPU
          amrex::Error("EOS: eos_input_tp has not been implemented gamma_law_2T EOS.");
#endif

          break;

      case eos_input_rp:

          // dens, pres, and xmass are inputs

          if constexpr (has_pressure<T>::value) {
              // Find the fraction of internal energy stored in heavies
              amrex::Real f = state.aux[0]/state.rho;
              // Calculate the electron temperature and set state.T to this value
              state.T = state.p * state.mu * m_nucleon * (1.0 - f) / (C::k_B * state.rho * x_e);
          }

          break;

      case eos_input_re:

          // dens, energy, and xmass are inputs
          // we already have the 2 temperatures from above

          if constexpr (has_energy<T>::value) {
              // Find the fraction of internal energy stored in heavies
              amrex::Real f = state.aux[0]/state.rho;
              // Calculate the electron temperature and set state.T to this value
              state.T = state.e * state.mu * m_nucleon * (1.0 - f) * (eos_rp::eos_gamma - 1.0) / (C::k_B * x_e);
          }

          break;

      case eos_input_ps:

          // pressure, entropy, and xmass are inputs

#ifndef AMREX_USE_GPU
          amrex::Error("EOS: eos_input_ps has not been implemented gamma_law_2T EOS.");
#endif

          break;

      case eos_input_ph:

          // pressure, enthalpy and xmass are inputs

#ifndef AMREX_USE_GPU
          amrex::Error("EOS: eos_input_ph has not been implemented gamma_law_2T EOS.");
#endif

          break;

      case eos_input_th:

          // temperature, enthalpy and xmass are inputs

          // This system is underconstrained.

  #ifndef AMREX_USE_GPU
          amrex::Error("EOS: eos_input_th is not a valid input for the gamma law EOS.");
  #endif

          break;

      default:

  #ifndef AMREX_USE_GPU
          amrex::Error("EOS: invalid input.");
  #endif

        break;

    }

    // Now we have the density and temperature (and mass fractions /
    // mu), regardless of the inputs.
    amrex::Real rhoinv = 1.0 / state.rho;

    // Compute the pressure simply from the ideal gas law, and the
    // specific internal energy using the gamma-law EOS relation.
    amrex::Real f = state.aux[0]/state.rho;
    amrex::Real pressure = state.rho * x_e * state.T * C::k_B / (state.mu * m_nucleon * (1.0 - f));
    amrex::Real energy = pressure / (eos_rp::eos_gamma - 1.0) * rhoinv;
    if constexpr (has_pressure<T>::value) {
        state.p = pressure;
    }
    if constexpr (has_energy<T>::value) {
        state.e = energy;
    }

    // enthalpy is h = e + p/rho
    if constexpr (has_enthalpy<T>::value) {
        state.h = energy + pressure * rhoinv;
    }

    // Compute the thermodynamic derivatives and specific heats
    amrex::Real Tinv = 1.0 / state.T;
    if constexpr (has_pressure<T>::value) {
        state.dpdT = state.p * Tinv;
        state.dpdr = state.p * rhoinv;
    }
    if constexpr (has_energy<T>::value) {
        state.dedT = state.e * Tinv;
        state.dedr = 0.0;
    }

    if constexpr (has_energy<T>::value) {
        state.cv = state.dedT;

        if constexpr (has_pressure<T>::value) {
            state.cp = eos_rp::eos_gamma * state.cv;

            state.gam1 = eos_rp::eos_gamma;

            state.dpdr_e = state.dpdr - state.dpdT * state.dedr * (1.0 / state.dedT);
            state.dpde = state.dpdT * (1.0 / state.dedT);

            // sound speed
            state.cs = std::sqrt(eos_rp::eos_gamma * state.p * rhoinv);
            if constexpr (has_G<T>::value) {
                state.G = 0.5 * (1.0 + eos_rp::eos_gamma);
            }
        }
    }

}


inline
void actual_eos_finalize() {

}

#endif
