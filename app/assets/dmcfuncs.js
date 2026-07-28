// Custom dmc.Tree node renderer for the filesystem browser.
//
// Registered on window.dashMantineFunctions and referenced from Python as
// renderNode={"function": "fsNode"}. Note that a custom renderNode makes
// Mantine ignore `checkboxes` / `expandedIcon` / `iconSide` — that is fine
// here, because only the *gating* tree needs checkboxes and it deliberately
// uses the built-in renderer.
//
// dash-iconify is used when present but is not required: the app must still
// render sensible icons if that bundle has not been loaded on the page.
var dmcfuncs = (window.dashMantineFunctions = window.dashMantineFunctions || {});

dmcfuncs.fsNode = function (payload) {
    var dmc = window.dash_mantine_components;
    var iconify = window.dash_iconify;
    var node = payload.node;
    var expanded = payload.expanded;
    var hasChildren = payload.hasChildren;
    var elementProps = payload.elementProps;

    var value = String(node.value || '');
    // Every real node value is an absolute path; fs_tree gives placeholder,
    // error and empty rows a "::prefix::" value instead. Those are shown dimmed
    // and italic so it is obvious they cannot be selected.
    var isSentinel = value.charAt(0) !== '/';

    function icon(name, color) {
        if (iconify && iconify.DashIconify) {
            return React.createElement(iconify.DashIconify, {
                icon: name,
                width: 16,
                height: 16,
                color: color
            });
        }
        // Fallback when the dash-iconify bundle is not on the page.
        var glyphs = {
            'tabler:folder': '📁',
            'tabler:folder-open': '📂',
            'tabler:file-analytics': '📄',
            'tabler:file': '📄',
            'tabler:dots': '…'
        };
        return React.createElement(
            'span',
            { style: { fontSize: 14 } },
            glyphs[name] || '•'
        );
    }

    function iconFor() {
        if (isSentinel) return icon('tabler:dots');
        if (hasChildren) {
            return icon(expanded ? 'tabler:folder-open' : 'tabler:folder', '#f59f00');
        }
        var lower = value.toLowerCase();
        if (lower.endsWith('.wsp')) return icon('tabler:file-analytics', '#4dabf7');
        if (lower.endsWith('.fcs')) return icon('tabler:file-analytics', '#69db7c');
        return icon('tabler:file', '#adb5bd');
    }

    var labelStyle = isSentinel
        ? { fontStyle: 'italic', opacity: 0.6 }
        : { whiteSpace: 'nowrap' };

    return React.createElement(
        dmc.Group,
        Object.assign({ key: value, gap: 6, wrap: 'nowrap' }, elementProps),
        iconFor(),
        React.createElement('span', { style: labelStyle }, node.label)
    );
};
