// Home page conda install box: clicking the command copies it (⧉ -> ✓).
// pkgdown picks this file up automatically, alongside pkgdown/extra.css.
document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('.ms-install .cmd').forEach(function (btn) {
    btn.addEventListener('click', function () {
      var cmd = btn.getAttribute('data-cmd');
      var flash = function () {
        btn.classList.add('done');
        setTimeout(function () { btn.classList.remove('done'); }, 1400);
      };
      // navigator.clipboard needs a secure context; fall back for file:// and http://
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(cmd).then(flash);
      } else {
        var ta = document.createElement('textarea');
        ta.value = cmd;
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        try { document.execCommand('copy'); flash(); } catch (e) { /* no-op */ }
        document.body.removeChild(ta);
      }
    });
  });
});
