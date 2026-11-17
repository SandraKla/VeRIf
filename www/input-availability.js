(function () {
  'use strict';

  var dataAvailable = false;
  var downloadsAvailable = false;
  var downloadIds = ['download_ritable', 'download_zlogtable'];
  var downloadUrls = {};

  function isAnalysisTab(link) {
    return link.matches('#analysis_tab > li > a') &&
      ['Input', 'About'].indexOf(link.getAttribute('data-value')) === -1;
  }

  function isLocked(link) {
    return (!dataAvailable && (isAnalysisTab(link) || link.matches('#results_tab > li > a'))) ||
      (!downloadsAvailable && downloadIds.indexOf(link.id) !== -1);
  }

  function syncDownloads() {
    downloadIds.forEach(function (id) {
      var link = document.getElementById(id);
      if (!link) return;
      var url = link.getAttribute('href');
      if (url && url !== '#') downloadUrls[id] = url;
      if (!downloadsAvailable) {
        if (url !== null) link.removeAttribute('href');
      } else if (!url && downloadUrls[id]) {
        link.setAttribute('href', downloadUrls[id]);
      }
      link.classList.toggle('disabled', !downloadsAvailable);
      link.setAttribute('aria-disabled', String(!downloadsAvailable));
      link.setAttribute('tabindex', downloadsAvailable ? '0' : '-1');
      link.title = downloadsAvailable ? '' : 'Available after a valid reference interval has been calculated.';
    });
  }

  function syncAvailability() {
    document.body.classList.toggle('verif-data-ready', dataAvailable);
    document.querySelectorAll('#analysis_tab > li > a, #results_tab > li > a').forEach(function (link) {
      var locked = isLocked(link);
      link.parentElement.classList.toggle('disabled', locked);
      link.setAttribute('aria-disabled', String(locked));
      if (locked) {
        link.setAttribute('tabindex', '-1');
        link.title = 'Load data in the Input tab first.';
      } else {
        link.removeAttribute('tabindex');
        link.removeAttribute('title');
      }
    });
    syncDownloads();
    var hint = document.getElementById('downloads_unavailable');
    if (hint) hint.hidden = downloadsAvailable;
    if (!dataAvailable) {
      var active = document.querySelector('#analysis_tab > li.active > a');
      if (active && isAnalysisTab(active)) {
        window.jQuery('#analysis_tab > li > a[data-value="Input"]').tab('show');
      }
    }
  }

  // Capture events before Bootstrap or Shiny can follow a disabled link.
  ['click', 'auxclick', 'keydown'].forEach(function (eventName) {
    document.addEventListener(eventName, function (event) {
      var link = event.target.closest('a');
      if (link && isLocked(link) &&
          (eventName !== 'keydown' || event.key === 'Enter' || event.key === ' ')) {
        event.preventDefault();
        event.stopImmediatePropagation();
      }
    }, true);
  });

  window.jQuery(function () {
    syncAvailability();
    window.jQuery('#analysis_tab, #results_tab').on('show.bs.tab', 'a', function (event) {
      if (isLocked(this)) event.preventDefault();
    });
    // Shiny supplies download URLs asynchronously; keep them inert until ready.
    var observer = new MutationObserver(syncDownloads);
    downloadIds.forEach(function (id) {
      var link = document.getElementById(id);
      if (link) observer.observe(link, { attributes: true, attributeFilter: ['href'] });
    });
  });

  Shiny.addCustomMessageHandler('setDataAvailability', function (message) {
    dataAvailable = message.dataAvailable === true;
    downloadsAvailable = dataAvailable && message.downloadsAvailable === true;
    syncAvailability();
  });
}());
