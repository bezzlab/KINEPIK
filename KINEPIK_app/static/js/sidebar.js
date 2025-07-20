document.addEventListener("DOMContentLoaded", function(){
    document.querySelectorAll(".section-toggle".forEach(button => {
        button.addEventListener("click",() => {
            const subsection = button.nextElementSibling;
            subsection.style.display = subsection.style.display === "block" ? "none" : "block";
        });
    }));

    document.querySelectorAll('.subsection a').forEach(link => {
    link.addEventListener('click', (e) => {
      e.preventDefault();
      const page = link.getAttribute('data-page');
      loadDocPage(page);
    });
  });

  function loadDocPage(page) {
    const doc = {
      'auth-login': '<h2>Login</h2><p>To login, use POST /login with email and password.</p>',
      'auth-logout': '<h2>Logout</h2><p>Send a DELETE request to /logout with your token.</p>',
      'get-users': '<h2>GET /users</h2><p>Returns a list of all users.</p>',
      'post-login': '<h2>POST /login</h2><p>Authenticate user credentials and return a token.</p>'
    };

    document.getElementById('doc-content').innerHTML = doc[page] || '<p>Page not found.</p>';
  }
});


