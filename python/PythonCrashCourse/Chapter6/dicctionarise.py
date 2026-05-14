#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Play with dictionaries
# =============================================================================


def aliens_info():
    alien0 = {"color": "green", "points": 5}
    print(alien0)
    print(alien0["color"])
    print(alien0["points"])
    alien0["x_position"] = 0
    alien0["y_position"] = 25
    print(alien0)
    alien0["girth"] = alien0.get("girth", 0)
    print(alien0)
    for key, value in alien0.items():
        print(f"\nKey: {key}")
        print(f"Value: {value}")


def rivers_info():
    river_dictionary = {
        "Amazonas": "Brazil",
        "Nilo": "Egipto",
        "Ganges": "India",
        "ThreeGorgesRiver": "China",
        "Duero": "España",
        "Ebro": "España",
        "Po": "Francia",
    }
    for key, value in river_dictionary.items():
        print(f"Did yoy know that the river {key} flows through {value} country?")
    print("\nNow lets review the rivers:")
    for name in river_dictionary.keys():
        print(name)
    print("\nNow lets review the countries:")
    for value in set(river_dictionary.values()):
        print(value)
    famous_countries = ["Brazil", "España", "Alemania", "Francia", "Turquía"]
    print()
    for country in famous_countries:
        if country not in set(river_dictionary.values()):
            print(f"You need to write a river for {country}")


def multiple_aliens():
    aliens = []
    for alien_number in range(30):
        new_alien = {"color": "green", "points": 5, "speed": "slow"}
        aliens.append(new_alien)
    for alien in aliens[:5]:
        print(alien)
    print("...")

    for alien in aliens[:3]:
        if alien["color"] == "green":
            alien["color"] = "yellow"
            alien["points"] = 10
            alien["speed"] = "medium"
        elif alien["color"] == "yellow":
            alien["color"] = "red"
            alien["points"] = 15
            alien["speed"] = "fast"

    for alien in aliens[:5]:
        print(alien)
    print("...")
    print(aliens[0]["color"].title())


def main() -> None:
    multiple_aliens()
    pass


if __name__ == "__main__":
    main()
