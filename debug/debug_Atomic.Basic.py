
if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src.Atomic import BasicP

    wave = 6.563E-5 # [cm]
    freq = BasicP.wave_to_freq( wave )
    assert wave == BasicP.freq_to_wave( freq )

    v = 100 * 1.0E+5 # 100km/s
    print(wave)
    Ds = BasicP.Dvlocity_to_Dshift(wave, v)
    print(f"Doppler shift of {wave*1.E+8:1.1f}[A] in {v*1E-5:1.1f}[km/s] is {Ds*1.E+8:1.1f}[A]")

    Te = 1.E+4 # [K]
    Vt = 5 * 1.0E+5 # 5km/s
    am = 1 # hydrogen
    Dw = BasicP.get_Doppler_width(wave, Te, Vt, am)
    print(f"Doppler width of {wave*1.E+8:1.1f}[A] of hydrogen with {Vt*1E-5:1.1f}[km/s] turbulent velocity : {Dw*1.E+8:1.1f}[A]")
