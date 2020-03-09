// XXX drop init
Jupyter.notebook.kernel.comm_manager.register_target('xena_comm',
    function(comm, msg) {
        // comm is the frontend comm instance
        // Register handlers for later messages:
        comm.on_msg(function(msg) {
            console.log('client got msg', msg);
            window.xenabrowser.window.postMessage(msg.content.data, "*"); // fix *
        });
        comm.on_close(function(msg) {console.log('client comm closed');});
        window.addEventListener("message", function(evt) {
            if (evt.origin === window.xenabrowser.url) {
                console.log("notebook window received", evt.data);
                comm.send(evt.data)
            }
        }, false);
    });
// note that sessionStorage for the domain is inherited from the notebook tab, which
// makes the state of the xena client a bit sticky.
window.xenabrowser.window = window.open(window.xenabrowser.url + '/')
