#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 30. Aug 2026
#
# Purpose:      Explores the web
# =============================================================================

INTERFACE = False

from selenium import webdriver

options = webdriver.ChromeOptions()
options.add_experimental_option("excludeSwitches", ["enable-automation"])
options.add_experimental_option("useAutomationExtension", False)
options.add_argument("--disable-blink-features=AutomationControlled")
options.add_argument("--disable-extensions")
# options.add_argument("--disable-gpu")

# without Interface
if INTERFACE:
    options.add_argument("--headless")

driver = webdriver.Chrome(options=options)


# driver.get("https://duckduckgo.com")
driver.get("https://bot.sannysoft.com")

print("navigator.webdriver:", driver.execute_script("return navigator.webdriver"))
print("userAgent:", driver.execute_script("return navigator.userAgent"))
print("platform:", driver.execute_script("return navigator.platform"))
print("languages:", driver.execute_script("return navigator.languages"))
print("plugins:", driver.execute_script("return navigator.plugins.length"))

input("Press Enter to closethe browser ... ")

driver.quit()
