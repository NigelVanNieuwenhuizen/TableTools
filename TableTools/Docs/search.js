
function filterSidebar() {
    const searchInput = document.getElementById('searchInput');
    const filter = searchInput.value.toUpperCase();
    const sidebar = document.querySelector('.sidebar');
    const links = sidebar.querySelectorAll('a');
    
    // 1. Filter all links
    links.forEach(link => {
        // Check only anchors with an href
        if (link.href) {
            const text = link.textContent.toUpperCase();
            if (text.indexOf(filter) > -1) {
                link.style.display = "";
            } else {
                link.style.display = "none";
            }
        }
    });

    // 2. Handle nested toolboxes (<details> elements)
    const toolboxes = sidebar.querySelectorAll('.toolbox-details');
    toolboxes.forEach(toolbox => {
        // Get all visible links inside the details block (excluding the summary itself)
        const allLinksInside = toolbox.querySelectorAll('a:not(.toolbox-summary)');
        let anyLinkMatches = false;
        
        // Check if any child link is visible
        allLinksInside.forEach(link => {
            if (link.style.display !== "none") {
                anyLinkMatches = true;
            }
        });

        // Get the summary element (the clickable title for the toolbox)
        const summary = toolbox.querySelector('.toolbox-summary');
        
        if (filter === '') {
            // If search is empty, collapse all toolboxes/objects by default.
            toolbox.open = false; 
            if (summary) summary.style.display = "";
        } else if (anyLinkMatches) {
            // If a match is found inside, expand the toolbox and show the summary title.
            toolbox.open = true;
            if (summary) summary.style.display = "";
        } else {
            // If no match is found inside, collapse it. 
            // Check if the summary title itself matches the filter.
            if (summary && summary.textContent.toUpperCase().includes(filter)) {
                 summary.style.display = "";
                 toolbox.open = true;
            } else {
                 toolbox.open = false;
                 if (summary) summary.style.display = "none";
            }
        }
    });
}
