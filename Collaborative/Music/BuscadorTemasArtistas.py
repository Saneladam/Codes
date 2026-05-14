import yt_dlp
import os

# --- CONFIGURACIÓN ---
ARCHIVO_ARTISTAS = 'artistas.txt'
CARPETA_DESTINO = 'Temas mp3'

def obtener_config_descarga(artista_actual):
    # Creamos la ruta completa: Temas mp3/Artista - Titulo.mp3
    ruta_salida = os.path.join(CARPETA_DESTINO, f'{artista_actual} - %(title)s.%(ext)s')
    
    return {
        'format': 'bestaudio/best',
        'postprocessors': [{
            'key': 'FFmpegExtractAudio',
            'preferredcodec': 'mp3',
            'preferredquality': '192',
        }],
        'outtmpl': ruta_salida,
        'quiet': False,
        'no_warnings': True,
    }

def ejecutar_proyecto():
    # 1. Verificar archivo de artistas
    if not os.path.exists(ARCHIVO_ARTISTAS):
        print(f"❌ Error: No encuentro el archivo '{ARCHIVO_ARTISTAS}'")
        return

    # 2. Crear la subcarpeta si no existe
    if not os.path.exists(CARPETA_DESTINO):
        os.makedirs(CARPETA_DESTINO)
        print(f"📁 Carpeta '{CARPETA_DESTINO}' creada con éxito.")

    with open(ARCHIVO_ARTISTAS, 'r', encoding='utf-8') as f:
        artistas = [linea.strip() for linea in f if linea.strip()]

    biblioteca = {}
    
    print("🔍 FASE 1: Analizando los 5 temas más virales por artista...")
    
    with yt_dlp.YoutubeDL({'quiet': True, 'extract_flat': True}) as ydl:
        for artista in artistas:
            print(f"  -> Analizando hits de: {artista}")
            busqueda = f"ytsearch15:{artista} official audio"
            try:
                info = ydl.extract_info(busqueda, download=False)
                videos = info.get('entries', [])
                
                candidatos = []
                for v in videos:
                    vistas = v.get('view_count')
                    if vistas:
                        candidatos.append({
                            'url': v.get('url'), 
                            'titulo': v.get('title'), 
                            'vistas': vistas
                        })
                
                top_3 = sorted(candidatos, key=lambda x: x['vistas'], reverse=True)[:10]
                biblioteca[artista] = top_3
                
            except Exception as e:
                print(f"     ⚠️ No se pudo analizar a {artista}: {e}")

    # --- FASE DE DESCARGA ---
    print("\n" + "="*50)
    print(f"RESUMEN: Se descargarán en la carpeta '{CARPETA_DESTINO}'")
    
    confirmar = input("\n¿Proceder con la descarga? (s/n): ").lower()

    if confirmar == 's':
        print("\n🚀 FASE 2: Iniciando conversión a MP3...")
        for artista, temas in biblioteca.items():
            with yt_dlp.YoutubeDL(obtener_config_descarga(artista)) as ydl_down:
                for tema in temas:
                    try:
                        print(f"\nBajando: {artista} - {tema['titulo']}")
                        ydl_down.download([tema['url']])
                    except Exception as e:
                        print(f"❌ Falló la descarga de {tema['titulo']}: {e}")
        
        print(f"\n✨ ¡Todo listo! Revisa la carpeta '{CARPETA_DESTINO}'.")
    else:
        print("Operación cancelada.")

if __name__ == "__main__":
    ejecutar_proyecto()