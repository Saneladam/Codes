#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 30. Aug 2026
#
# Purpose:      Explores the web
# =============================================================================

from selenium import webdriver

options = webdriver.ChromeOptions()
# options.add_experimental_option("exludeSwitches", ["enable-automation"])
# options.add_experimental_option("useAutomationExtension", False)
options.add_argument("--disable-blinnk-features=AutomationControlled")

driver = webdriver.Chrome()


# driver.get("https://duckduckgo.com")
driver.get("https://bot.sannysoft.com")

input("Press Enter to closethe browser ... ")

driver.quit()
