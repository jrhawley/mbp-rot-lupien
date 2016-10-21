module.exports = function(grunt) {

    // Project configuration.
    grunt.initConfig({
        pkg: grunt.file.readJSON('package.json'),
        assemble: {
            options: {
                assets: "path/to/assets",
                data:   "path/to/config.json" 
            },
            project: {
                options: {
                    layout: "path/to/default-layout.hbs",
                    partials: "path/to/partials/**/*.hbs" 
                },
                files: {
                    'dest': ["path/to/pages/**/*.hbs" ]
                }
            }
        }
    });

    // Load assemble markdown task
    grunt.loadNpmTasks('grunt-assemble');

    // Default task(s).
    grunt.registerTask('default', ['assemble']);

};