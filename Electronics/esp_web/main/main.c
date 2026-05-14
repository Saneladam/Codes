#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "esp_wifi.h"
#include "esp_event.h"
#include "nvs_flash.h"
#include "esp_http_server.h"

static const char* html_page =
"<!DOCTYPE html>"
"<html><body>"
"<h1>ESP32-C3 Minimal</h1>"
"</body></html>";

esp_err_t root_handler(httpd_req_t *req)
{
    httpd_resp_send(req, html_page, HTTPD_RESP_USE_STRLEN);
    return ESP_OK;
}

void app_main(void)
{
    // Inicializar NVS
    nvs_flash_init();

    // Aquí iría init WiFi (luego lo hacemos)
}
