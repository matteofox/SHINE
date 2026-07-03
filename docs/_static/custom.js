document.addEventListener("DOMContentLoaded", function() {
    // Create the Lightbox elements dynamically
    const lightboxModal = document.createElement("div");
    lightboxModal.className = "shine-lightbox-modal";
    lightboxModal.id = "shineLightbox";
    
    const closeBtn = document.createElement("span");
    closeBtn.className = "shine-lightbox-close";
    closeBtn.innerHTML = "&times;";
    
    const lightboxImg = document.createElement("img");
    lightboxImg.className = "shine-lightbox-content";
    lightboxImg.id = "shineLightboxImg";
    
    lightboxModal.appendChild(closeBtn);
    lightboxModal.appendChild(lightboxImg);
    document.body.appendChild(lightboxModal);
    
    // Attach click events to images in the main content body
    const contentImages = document.querySelectorAll(".rst-content img");
    contentImages.forEach(img => {
        // Skip small icon images if any
        if (img.width < 50 || img.height < 50) return;
        
        img.classList.add("zoomable-img");
        img.addEventListener("click", function(e) {
            e.preventDefault();
            lightboxImg.src = this.src;
            lightboxModal.style.display = "flex";
            document.body.style.overflow = "hidden"; // Disable background scrolling
        });
    });
    
    // Function to close the Lightbox
    function closeLightbox() {
        lightboxModal.style.display = "none";
        document.body.style.overflow = "auto"; // Re-enable background scrolling
    }
    
    // Close events
    lightboxModal.addEventListener("click", function(e) {
        if (e.target !== lightboxImg) {
            closeLightbox();
        }
    });
    
    closeBtn.addEventListener("click", closeLightbox);
    
    document.addEventListener("keydown", function(e) {
        if (e.key === "Escape" && lightboxModal.style.display === "flex") {
            closeLightbox();
        }
    });
});
