// Load available tools on page load
document.addEventListener('DOMContentLoaded', async () => {
    await loadTools();
    await checkHealth();
    
    // Enable Enter key to send message
    document.getElementById('user-input').addEventListener('keypress', (e) => {
        if (e.key === 'Enter' && !e.shiftKey) {
            e.preventDefault();
            sendMessage();
        }
    });
});

async function checkHealth() {
    try {
        const response = await fetch('/health');
        const health = await response.json();
        console.log('System health:', health);
    } catch (error) {
        console.error('Health check failed:', error);
        showError('System health check failed. Please check if all services are running.');
    }
}

async function loadTools() {
    try {
        const response = await fetch('/tools');
        const tools = await response.json();
        
        const toolsList = document.getElementById('tools-list');
        toolsList.innerHTML = '';
        
        if (tools.length === 0) {
            toolsList.innerHTML = '<p class="no-tools">No tools available. Check BioMCP connection.</p>';
            return;
        }
        
        tools.forEach(tool => {
            const toolItem = document.createElement('div');
            toolItem.className = 'tool-item';
            
            let parametersHtml = '';
            if (tool.parameters && tool.parameters.length > 0) {
                const params = tool.parameters.map(p => 
                    `<span class="param">${p.name}${p.required ? '*' : ''}</span>`
                ).join(', ');
                parametersHtml = `<div class="tool-params">Parameters: ${params}</div>`;
            }
            
            toolItem.innerHTML = `
                <div class="tool-name">${tool.name}</div>
                <div class="tool-description">${tool.description || 'No description available'}</div>
                ${parametersHtml}
            `;
            toolsList.appendChild(toolItem);
        });
    } catch (error) {
        console.error('Error loading tools:', error);
        document.getElementById('tools-list').innerHTML = '<p class="error">Error loading tools</p>';
    }
}

async function sendMessage() {
    const input = document.getElementById('user-input');
    const message = input.value.trim();
    
    if (!message) return;
    
    // Disable input while processing
    input.disabled = true;
    document.getElementById('send-btn').disabled = true;
    
    // Add user message to chat
    addMessage(message, 'user');
    
    // Clear input
    input.value = '';
    
    // Add loading indicator
    const loadingId = addMessage('<div class="loading"></div>', 'assistant');
    
    try {
        const response = await fetch('/chat', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ message: message })
        });
        
        const reader = response.body.getReader();
        const decoder = new TextDecoder();
        
        let assistantMessage = '';
        const messageElement = document.getElementById(loadingId);
        let isFirstChunk = true;
        
        while (true) {
            const { done, value } = await reader.read();
            if (done) break;
            
            const chunk = decoder.decode(value);
            const lines = chunk.split('\n');
            
            for (const line of lines) {
                if (line.startsWith('data: ')) {
                    try {
                        const data = JSON.parse(line.slice(6));
                        
                        if (isFirstChunk) {
                            messageElement.innerHTML = '';
                            isFirstChunk = false;
                        }
                        
                        switch (data.type) {
                            case 'content':
                                assistantMessage += data.content;
                                messageElement.innerHTML = formatMessage(assistantMessage);
                                break;
                                
                            case 'tool_call':
                                addToolCall(data.tool, data.parameters);
                                break;
                                
                            case 'tool_result':
                                addToolResult(data.tool, data.result);
                                break;
                                
                            case 'interpretation':
                                assistantMessage += data.content;
                                messageElement.innerHTML = formatMessage(assistantMessage);
                                break;
                                
                            case 'tool_error':
                                addToolError(data.tool, data.error);
                                break;
                        }
                        
                        // Auto-scroll to bottom
                        const chatMessages = document.getElementById('chat-messages');
                        chatMessages.scrollTop = chatMessages.scrollHeight;
                        
                    } catch (e) {
                        console.error('Error parsing chunk:', e);
                    }
                }
            }
        }
    } catch (error) {
        console.error('Error:', error);
        document.getElementById(loadingId).innerHTML = 'Error: Failed to get response';
    } finally {
        // Re-enable input
        input.disabled = false;
        document.getElementById('send-btn').disabled = false;
        input.focus();
    }
}

function addMessage(content, sender) {
    const messagesContainer = document.getElementById('chat-messages');
    const messageDiv = document.createElement('div');
    const messageId = `msg-${Date.now()}`;
    
    messageDiv.id = messageId;
    messageDiv.className = `message ${sender}-message`;
    messageDiv.innerHTML = sender === 'user' ? content : formatMessage(content);
    
    messagesContainer.appendChild(messageDiv);
    messagesContainer.scrollTop = messagesContainer.scrollHeight;
    
    return messageId;
}

function addToolCall(toolName, parameters) {
    const messagesContainer = document.getElementById('chat-messages');
    const toolDiv = document.createElement('div');
    toolDiv.className = 'tool-call';
    toolDiv.innerHTML = `
        <div class="tool-header">🔧 Calling Tool: ${toolName}</div>
        <div class="tool-params-display">${JSON.stringify(parameters, null, 2)}</div>
    `;
    messagesContainer.appendChild(toolDiv);
}

function addToolResult(toolName, result) {
    const messagesContainer = document.getElementById('chat-messages');
    const resultDiv = document.createElement('div');
    resultDiv.className = 'tool-result';
    resultDiv.innerHTML = `
        <div class="tool-header">✅ Tool Result: ${toolName}</div>
        <div class="tool-result-content">${formatMessage(result)}</div>
    `;
    messagesContainer.appendChild(resultDiv);
}

function addToolError(toolName, error) {
    const messagesContainer = document.getElementById('chat-messages');
    const errorDiv = document.createElement('div');
    errorDiv.className = 'tool-error';
    errorDiv.innerHTML = `
        <div class="tool-header">❌ Tool Error: ${toolName}</div>
        <div class="tool-error-content">${error}</div>
    `;
    messagesContainer.appendChild(errorDiv);
}

function formatMessage(content) {
    // Convert markdown-like formatting to HTML
    return content
        .replace(/```([\s\S]*?)```/g, '<pre><code>$1</code></pre>')
        .replace(/`([^`]+)`/g, '<code>$1</code>')
        .replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>')
        .replace(/\n/g, '<br>');
}

function showError(message) {
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error-notification';
    errorDiv.textContent = message;
    document.body.appendChild(errorDiv);
    
    setTimeout(() => {
        errorDiv.remove();
    }, 5000);
}

