/*
 * Click-to-enlarge lightbox, shared across eartigau's tools.
 * Pairs with the #lightbox markup in banner-snippet.html and the
 * .lightbox* rules in theme.css.
 *
 * Works via event delegation on `document`, so it needs no per-image
 * wiring and picks up images added to the DOM later (e.g. a gallery
 * built from a fetched manifest). Caption resolution order:
 *   1. the clicked image's own <figure><figcaption>, if any
 *   2. its enclosing .card's <h2> + .hint text
 *   3. the image's alt text
 *
 * Call setupLightbox() once after the #lightbox markup exists in the DOM
 * (a plain <script> at the end of <body> is fine, no need to wait for
 * DOMContentLoaded).
 */

function captionFor(img) {
  const figure = img.closest("figure");
  if (figure) {
    const cap = figure.querySelector("figcaption");
    if (cap && cap.textContent.trim()) {
      return cap.textContent.trim();
    }
  }
  const card = img.closest(".card");
  if (card) {
    const h2 = card.querySelector("h2");
    const hint = card.querySelector(".hint");
    const title = h2 ? h2.textContent.trim() : "";
    const desc = hint ? hint.textContent.trim() : "";
    return desc ? `${title}: ${desc}` : title;
  }
  return img.alt || "";
}

function setupLightbox() {
  const lightbox = document.getElementById("lightbox");
  const lightboxImg = document.getElementById("lightboxImg");
  const lightboxCaption = document.getElementById("lightboxCaption");
  const lightboxClose = document.getElementById("lightboxClose");
  if (!lightbox || !lightboxImg || !lightboxCaption || !lightboxClose) {
    return;
  }

  function open(img) {
    lightboxImg.src = img.currentSrc || img.src;
    lightboxImg.alt = img.alt || "";
    lightboxCaption.textContent = captionFor(img);
    lightbox.hidden = false;
  }

  function close() {
    lightbox.hidden = true;
    lightboxImg.src = "";
  }

  document.addEventListener("click", (e) => {
    const img = e.target.closest(".image-frame img, .gallery-item img, .lecture-figure img");
    if (img) {
      open(img);
    }
  });

  lightbox.addEventListener("click", (e) => {
    if (e.target === lightbox) {
      close();
    }
  });
  lightboxClose.addEventListener("click", close);
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && !lightbox.hidden) {
      close();
    }
  });
}

setupLightbox();
