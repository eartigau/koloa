/*
 * koloa's site: one page at a time.
 *
 * The content is a list of pages (<section class="view" id=...>); the
 * sidebar links to them by their id (#id). This script shows the page the
 * address names (the overview when it names none), marks it in the sidebar,
 * opens its group there (use cases, modules, mathematics) and puts links to
 * the previous and next pages, in the order of the sidebar, under it. The
 * back and forward buttons of the browser move between pages. A link to an
 * element inside a page (#some-id) shows that page and scrolls to it.
 *
 * Without this script every page shows, one under another. It also picks
 * the photo of the overview at random.
 */
(function () {
  "use strict";

  /* the photo of the overview: one of the koloa photos, at random (without
     this script, the first) */
  var picks = document.querySelectorAll(".hero-pick");
  if (picks.length > 1) {
    var chosen = Math.floor(Math.random() * picks.length);
    Array.prototype.forEach.call(picks, function (pick, it) {
      pick.hidden = it !== chosen;
    });
  }

  var views = Array.prototype.slice.call(document.querySelectorAll(".view"));
  if (!views.length) {
    return;
  }
  var nav = document.querySelector(".sidebar");
  var details = nav ? nav.querySelector("details") : null;
  var pager = document.getElementById("pager");
  var narrow = window.matchMedia("(max-width: 900px)");
  var home = views[0].id;
  var baseTitle = document.title;

  /* the pages in the order of the sidebar (the order of the page) */
  var order = views.map(function (view) {
    return {id: view.id, title: view.getAttribute("data-title") || view.id};
  });

  function viewOf(id) {
    var target = id ? document.getElementById(id) : null;
    if (!target) {
      return null;
    }
    return target.classList.contains("view") ? target : target.closest(".view");
  }

  function markSidebar(view) {
    if (!nav) {
      return;
    }
    var group = view.getAttribute("data-group") || view.id;
    Array.prototype.forEach.call(nav.querySelectorAll("a[href^='#']"), function (link) {
      link.classList.toggle("active", link.getAttribute("href") === "#" + view.id);
    });
    Array.prototype.forEach.call(nav.querySelectorAll(".nav-group"), function (item) {
      item.classList.toggle("open", item.getAttribute("data-group") === group);
    });
    var active = nav.querySelector("a.active");
    if (active && !narrow.matches && active.scrollIntoView) {
      active.scrollIntoView({block: "nearest"});
    }
  }

  function fillPager(view) {
    if (!pager) {
      return;
    }
    var it = order.map(function (page) { return page.id; }).indexOf(view.id);
    var prev = it > 0 ? order[it - 1] : null;
    var next = it >= 0 && it < order.length - 1 ? order[it + 1] : null;
    pager.innerHTML = "";
    [[prev, "prev", "← "], [next, "next", ""]].forEach(function (spec) {
      var page = spec[0];
      var link = document.createElement("a");
      link.className = "pager-" + spec[1];
      if (page) {
        link.href = "#" + page.id;
        link.textContent = spec[1] === "prev" ? spec[2] + page.title : page.title + " →";
      } else {
        link.setAttribute("aria-hidden", "true");
      }
      pager.appendChild(link);
    });
  }

  function show(id) {
    var view = viewOf(id) || document.getElementById(home);
    views.forEach(function (page) {
      page.hidden = page !== view;
    });
    markSidebar(view);
    fillPager(view);
    /* the other language, on the same page */
    Array.prototype.forEach.call(document.querySelectorAll("a.lang-link"),
      function (link) {
        link.href = link.getAttribute("href").split("#")[0] + "#" + view.id;
      });
    var title = view.getAttribute("data-title");
    document.title = view.id === home || !title ? baseTitle : title + " · koloa";
    /* the top of the page, or the element the address names */
    var target = id ? document.getElementById(id) : null;
    if (target && target !== view) {
      target.scrollIntoView({block: "start"});
    } else {
      /* after the jump of the browser to the anchor, back to the top */
      window.scrollTo(0, 0);
      window.setTimeout(function () {
        window.scrollTo(0, 0);
      }, 0);
    }
  }

  document.documentElement.classList.add("paged");
  window.addEventListener("hashchange", function () {
    show(window.location.hash.slice(1));
  });
  show(window.location.hash.slice(1));

  /* the sidebar folds into a menu on a narrow screen, and closes after a
     choice there */
  function fold() {
    if (details) {
      details.open = !narrow.matches;
    }
  }
  fold();
  if (narrow.addEventListener) {
    narrow.addEventListener("change", fold);
  }
  if (nav) {
    nav.addEventListener("click", function (event) {
      if (event.target.closest && event.target.closest("a") && narrow.matches && details) {
        details.open = false;
      }
    });
  }
})();
