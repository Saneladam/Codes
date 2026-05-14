#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      if statements
# =============================================================================

bitches = ["Natalia", "Carlota", "Rebeca", "Alejandra", "Carla", "Sandra"]
sifilis = ["Marta", "Elena", "Laura", "Rut", "Adriana"]
# bitches = []


def main() -> None:
    if bitches:
        guarrilla = input("Cómo te llamas descoñocida?\n> ")
        while guarrilla in sifilis:
            print("No soy yo, eres tu... es que eres lesbiana")
            guarrilla = input("Cómo te llamas descoñocida?\n> ")
        if guarrilla not in bitches:
            print("Lo siento puta, no eres bien recibida")
        else:
            print("Eyyy que pasa zorri?")

        pass
    else:
        print("Tu prima esta buena")
    print("Recuerda, nunca te acerques a las gonorreicas de:\n")
    n = 0
    for sifilitica in sifilis:
        n += 1
        print(f"\t{n}:\t{sifilitica.upper()}")
    print("Y eso seria todo migo")


if __name__ == "__main__":
    main()
