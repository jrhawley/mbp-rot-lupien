$(document).ready(function() {
	/* Scroll-to-top button */
	// browser window scroll (in pixels) after which the "back to top" link is shown
	var offset = 300;
	//browser window scroll (in pixels) after which the "back to top" link opacity is reduced
	var offset_opacity = 1200;
	//duration of the top scrolling animation (in ms)
	var scroll_duration = 700;
	//grab the "back to top" and "back to bottom" links
	var back_to_top = $('.cd-top');
	var back_to_bottom = $('.cd-bottom');

	//hide or show the "back to top" link
	$(window).scroll(function() {
		var top_height_offset = $(this).scrollTop();
		(top_height_offset > offset) ? back_to_top.addClass('cd-is-visible'): back_to_top.removeClass('cd-is-visible');
	});

	//hide or show the "back to bottom" link
	$(window).scroll(function() {
		var bottom_height_offset = $(document).height()-$(window).height() - $(this).scrollTop();
		(bottom_height_offset > offset) ? back_to_bottom.addClass('cd-is-visible'): back_to_bottom.removeClass('cd-is-visible');
	});

	//smooth scroll to top
	back_to_top.on('click', function(event) {
		event.preventDefault();
		$('body,html').animate({
			scrollTop: 0,
		}, scroll_duration);
	});

	//smooth scroll to top
	back_to_bottom.on('click', function(event) {
		event.preventDefault();
		$('body,html').animate({
			scrollTop: $(document).height()-$(window).height(),
		}, scroll_duration);
	});

	$('.scrollspy').scrollSpy();
});