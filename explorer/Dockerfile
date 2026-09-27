# Exercise KG Explorer: static site (index.html + data.json) served by nginx.
# Everything is local: no CDN is needed at run time.
FROM nginx:1.27-alpine

COPY nginx.conf /etc/nginx/conf.d/default.conf
COPY vendor/ /usr/share/nginx/html/vendor/
COPY data.json /usr/share/nginx/html/data.json
COPY index.html /tmp/index.src.html

# index.html is written as a Claude artifact body (no <!doctype>/<meta charset>, cytoscape from a CDN).
# Turn it into a standalone page that loads the vendored cytoscape.
RUN { printf '<!doctype html>\n<html lang="en">\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n'; \
      sed 's#https://cdn.jsdelivr.net/npm/cytoscape@3.30.2/dist/cytoscape.min.js#vendor/cytoscape.min.js#' /tmp/index.src.html; \
      printf '\n</html>\n'; } > /usr/share/nginx/html/index.html \
 && grep -q 'vendor/cytoscape.min.js' /usr/share/nginx/html/index.html \
 && rm /tmp/index.src.html

EXPOSE 80
HEALTHCHECK --interval=30s --timeout=3s CMD wget -qO- http://127.0.0.1/healthz || exit 1
