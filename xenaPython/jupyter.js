Jupyter.notebook.kernel.comm_manager.register_target('xena_comm',
    (function() {
        var windows = [];
        return function(comm, msg) {
            // comm is the frontend comm instance
            // Register handlers for later messages:
            comm.on_msg(function(msg) {
                //console.log('client got msg', JSON.stringify(msg.content));
                var [w, id] = msg.content.data.id.split(/:/);
                msg.content.data.id = parseInt(id, 10);
                windows[w].postMessage(msg.content.data, "*"); // fix the "*"
            });
            comm.on_close(function(msg) {console.log('client comm closed');});
            window.addEventListener("message", function(evt) {
                if (evt.origin === window.xenabrowser.url) {
                    var w = windows.indexOf(evt.source);
                    if (w === -1) {
                        w = windows.length;
                        windows.push(evt.source);
                    }
                    //console.log("notebook window received", JSON.stringify(evt.data));
                    // prepend window to message id
                    evt.data.id = w.toString() + ':' + evt.data.id.toString();
                    comm.send(evt.data);
                }
            }, false);
        };
    }()));
// note that sessionStorage for the domain is inherited from the notebook tab, which
// makes the state of the xena client a bit sticky.
