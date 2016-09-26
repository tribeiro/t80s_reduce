
NBIAS = 11
NFLATS = 7
FILTERS = ['F378', 'G', 'I', 'F430', 'F861', 'F515', 'R', 'U', 'F660', 'Z', 'F410', 'F395']
DOWNLOAD_STRING_HEADER = 'wget -c --no-check-cert '
DOWNLOAD_STRING_HOST = 'https://t80s_images:t80s_images_keywords_pass@splus.astro.ufsc.br//{night}/{filename}' \
                       ' -P {dest_path}'